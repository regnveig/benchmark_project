__author__ = "Polina Belokopytova & Emil Viesná"
__version__ = "0.1b"
__date__ = "2021 Nov 1"

from contextlib import contextmanager
from cooltools.insulation import calculate_insulation_score
from hicreppy import hicrep
from multiprocessing import cpu_count
from sklearn.metrics import precision_recall_curve, auc
import argparse
import cooler
import datetime
import hashlib
import json
import logging
import matplotlib.pyplot as plt
import numpy
import os, sys
import pandas, pandarallel
import subprocess
import tempfile
import time
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)
pandarallel.pandarallel.initialize(nb_workers = cpu_count(), verbose = 0)

# ------======| CONST |======------

C_SCC_MAXDIST = 1500000
C_SCC_H = 2
C_CONTACT_COEF = 100
C_RANDOM_INTER_N = 5000
C_RANDOM_INTER_SIGMA = 2

## ------======| LOGGING |======------

def ExceptionHook(Type, Value, Traceback): logging.exception(f"{Type.__name__}: {Value}", exc_info=(Type, Value, Traceback))

def ConfigureLogger(LogFileName: str = os.devnull, Level: int = logging.INFO) -> None:
	Formatter = "%(asctime)-30s%(levelname)-13s%(funcName)-35s%(message)s"
	logging.basicConfig(level = Level, format = Formatter, handlers = [logging.FileHandler(LogFileName), logging.StreamHandler(sys.stderr)], force = True)
	sys.excepthook = ExceptionHook

def SecToTime(Sec: float) -> str: return f"{int(Sec / 3600):02}:{int((Sec // 60) % 60):02}:{int(Sec % 60):02}"

@contextmanager
def Timer(EndMessage: str, StartMessage: str = None):
	if StartMessage is not None: logging.info(StartMessage)
	StartTime = time.time()
	yield
	logging.info(f"{EndMessage} [{SecToTime(time.time() - StartTime)}]")

# ------======| DATA LOADING |======------

def BinSearch(Chrom, End, BinDict):
	try:
		return BinDict[(Chrom, End)]
	except KeyError:
		raise RuntimeError(f"Unknown bin")

def Tsv2Cool(TsvFN, OutputCoolFN, TemplateCoolFN, Chrom, BinSize):
	for Line in [f"Input TSV: {TsvFN}", f"Output COOL: {OutputCoolFN}", f"Template COOL: {TemplateCoolFN}", f"Chrom: {Chrom}", f"Resolution: {int(BinSize / 1000)} kb"]: logging.info(Line)
	DType = {"chrom": str, "end1": int, "end2": int, "balanced": float}
	Pixels = pandas.read_csv(TsvFN, sep='\t', names=DType.keys(), dtype=DType, header=None)
	Bins = cooler.Cooler(TemplateCoolFN).bins().fetch(Chrom)
	Bins["weight"] = 1 / C_CONTACT_COEF
	BinsDict = {(Line["chrom"], Line["end"]): index for index, Line in Bins.iterrows()}
	for num in [1, 2]: Pixels[f"bin{num}_id"] = Pixels.parallel_apply(lambda x: BinSearch(x["chrom"], x[f"end{num}"], BinsDict), axis=1)
	Pixels["count"] = Pixels["balanced"].parallel_apply(lambda x: int(x * pow(C_CONTACT_COEF, 2)))
	Pixels = Pixels[["bin1_id", "bin2_id", "count"]]
	cooler.create_cooler(OutputCoolFN, Bins, Pixels)

def Cool2Cool(InputCoolFN, OutputCoolFN, Chrom):
	for Line in [f"Input COOL: {InputCoolFN}", f"Output COOL: {OutputCoolFN}", f"Chrom: {Chrom}"]: logging.info(Line)
	Data = cooler.Cooler(InputCoolFN)
	Bins = Data.bins().fetch(Chrom)
	Pixels = Data.matrix(as_pixels=True, balance=False).fetch(Chrom)
	cooler.create_cooler(OutputCoolFN, Bins, Pixels)

def AlignCools(InputFNA, InputFNB, OutputFNA, OutputFNB, Chrom):
	for Line in [f"Input COOL [A]: {InputFNA}", f"Output COOL [A]: {OutputFNA}", f"Input COOL [B]: {InputFNB}", f"Output COOL [B]: {OutputFNB}", f"Chrom: {Chrom}"]: logging.info(Line)
	InputA, InputB = cooler.Cooler(InputFNA), cooler.Cooler(InputFNB)
	BinsA, BinsB = InputA.bins().fetch(Chrom), InputB.bins().fetch(Chrom)
	PixelsA, PixelsB = InputA.matrix(as_pixels=True, balance=False).fetch(Chrom), InputB.matrix(as_pixels=True, balance=False).fetch(Chrom)
	MergePixels = pandas.merge(PixelsA, PixelsB, how="inner", on=["bin1_id", "bin2_id"])
	PixelsA, PixelsB = MergePixels[["bin1_id", "bin2_id", "count_x"]].rename(columns={"count_x": "count"}), MergePixels[["bin1_id", "bin2_id", "count_y"]].rename(columns={"count_y": "count"})
	cooler.create_cooler(OutputFNA, BinsA, PixelsA)
	cooler.create_cooler(OutputFNB, BinsB, PixelsB)

def GetMatrix(Cool, Chrom): return Cool.matrix(as_pixels=True, balance=True).fetch(Chrom)

def InsulationData(Datasets, Window):
	for Line in [f"Window: {Window}"]: logging.info(Line)
	InsScores = {Type: calculate_insulation_score(Data, [Window], ignore_diags=2, append_raw_scores=True).rename(columns={f"sum_balanced_{Window}": f"sum_balanced_{'-'.join(Type)}"}) for Type, Data in Datasets.items()}
	Result = None
	for Type, Data in InsScores.items(): Result = Data.copy() if Result is None else pandas.merge(Result, Data, how="inner", on=["chrom", "start", "end"])
	Result = Result[["chrom", "start", "end"] + [f"sum_balanced_{'-'.join(Type)}" for Type in InsScores.keys()]]
	for DT in ["Exp", "Pred"]: 
		Result[f"sum_balanced_Mut/Wt-{DT}"] = Result.apply(lambda x: 0 if x[f"sum_balanced_Wt-{DT}"] == 0 else (x[f"sum_balanced_Mut-{DT}"] / x[f"sum_balanced_Wt-{DT}"]), axis=1)
		NonZero = Result[f"sum_balanced_Mut/Wt-{DT}"][Result[f"sum_balanced_Mut/Wt-{DT}"] != 0]
		Mean, Std = numpy.mean(NonZero), numpy.std(NonZero)
		NonZero = NonZero.apply(lambda x: (x - Mean) / Std)
		NonZero.name = f"sum_balanced_sigma_Mut/Wt-{DT}"
		Result = pandas.concat([Result, NonZero], axis=1)
	Result["Y-True"] = Result["sum_balanced_sigma_Mut/Wt-Exp"].apply(lambda x: (x == x) and ((x > 2) or (x < -2)))
	return Result

def EctopicInteractionsArray(CoolWT, CoolMut, Chrom, CaptureStart, CaptureEnd, RearrStart, RearrEnd, Normalized):
	Capture = f"{Chrom}:{CaptureStart}-{CaptureEnd}"
	Rearr = f"{Chrom}:{RearrStart}-{RearrEnd}"
	for Line in [f"WT COOL: {CoolWT.store}", f"Mut COOL: {CoolMut.store}", f"Capture: {Capture}", f"Rearrangement: {Rearr}", f"Normalized: {'yes' if Normalized else 'no'}"]: logging.info(Line)
	
	def PrepareData(Cool):
		Data = numpy.nan_to_num(Cool.matrix(balance=False).fetch(Capture))
		RearrStartBin = (RearrStart - CaptureStart) // Cool.binsize
		RearrEndBin = (RearrEnd - CaptureStart) // Cool.binsize
		Data[RearrStartBin:RearrEndBin+1, :] = numpy.zeros(Data[RearrStartBin:RearrEndBin+1, :].shape)
		Data[:, RearrStartBin:RearrEndBin+1] = numpy.zeros(Data[:, RearrStartBin:RearrEndBin+1].shape)
		return Data
	
	SumData = lambda Data: sum(list(map(sum, Data)))
	
	DataWT, DataMut = PrepareData(CoolWT), PrepareData(CoolMut)
	SumWT, SumMut = SumData(DataWT), SumData(DataMut)
	DataMut = DataMut * (1 if Normalized else (SumWT / SumMut))
	DiffArray = DataMut - DataWT
	
	logging.info(f"Diff array prepared")
	
	for i in range(len(DiffArray)):
		X, Y = numpy.array(range(0, len(DiffArray) - i)), numpy.array(range(i, len(DiffArray)))
		Coeff = numpy.average(DataWT[X,Y])
		if Coeff != 0: DiffArray[X,Y], DiffArray[Y,X] = (DiffArray[X,Y] / Coeff), (DiffArray[X,Y] / Coeff)
	
	logging.info(f"Diff matrix normalized")
	
	DiffDiags = [(k, numpy.diag(DiffArray, k=k)) for k in range(len(DiffArray))]
	numpy.nan_to_num(DiffDiags, copy=False)
	
	EctopicArray = numpy.zeros_like(DataMut)
	
	logging.info(f"Diagonals created")
	
	for k, kDiag in DiffDiags:
		NonZero = numpy.nonzero(kDiag)
		if NonZero[0].size == 0: continue
		PercTop, PercBottom = numpy.percentile(kDiag[NonZero], 96), numpy.percentile(kDiag[NonZero], 4)
		Diag96Perc = [item for item in kDiag if (PercBottom < item < PercTop) and (item != 0)]
		if len(Diag96Perc) < 10: continue
		DiagMean, DiagStd = numpy.mean(Diag96Perc), numpy.std(Diag96Perc)
		if DiagStd == 0: continue
		for n, Contact in enumerate(kDiag):
			if Contact != 0: EctopicArray[k + n, n] = (Contact - DiagMean) / DiagStd
	
	logging.info(f"Ectopic interactions on diagonals found")
	
	return EctopicArray

def IntersectEctopicMatrices(MatrixA, MatrixB, SD):
	Condition = lambda Matrix: numpy.logical_and(numpy.isfinite(Matrix), numpy.logical_or(Matrix < -SD, Matrix > SD))
	return numpy.sum(numpy.logical_and(Condition(MatrixA), Condition(MatrixB)))

def MakeMcool(InputCool, OutputMcool, Resolution):
	for Line in [f"Input COOL: {InputCool}", f"Output MCOOL: {OutputMcool}", f"Resolution: {int(Resolution / 1000)} kb"]: logging.info(Line)
	cooler.zoomify_cooler(InputCool, OutputMcool, [Resolution], 1000)

# ------======| METRICS |======------

def PearsonCorr(SeriesA, SeriesB): return SeriesA.corr(SeriesB, method="pearson")

def SCC(CoolA, CoolB, MaxDist, h): return hicrep.genome_scc(CoolA, CoolB, max_dist=MaxDist, h=h)

def PRCurve(YTrue, Probas):
	YTrueNumpy = numpy.nan_to_num(YTrue)
	ProbasNumpy = numpy.nan_to_num(Probas)
	Precision, Recall, Thresholds = precision_recall_curve(YTrueNumpy, ProbasNumpy)
	
	return {
		"AUC": auc(Recall, Precision),
		"Precision": json.dumps(Precision.tolist()),
		"Recall": json.dumps(Recall.tolist()),
		"Thresholds": json.dumps(Thresholds.tolist())
		}

def RandomEctopicIntersections(EctopicArrayExp, EctopicArrayPred, n = C_RANDOM_INTER_N, sigma = C_RANDOM_INTER_SIGMA):
	for Line in [f"N: {n}", f"Sigma: {sigma}"]: logging.info(Line)
	RandomIntersections = []
	EctopicPredRandom = numpy.copy(EctopicArrayPred)
	Finite = numpy.where(numpy.isfinite(EctopicArrayPred))
	Rand = EctopicArrayPred[Finite].flatten()
	for i in range(n):
		numpy.random.shuffle(Rand)
		EctopicPredRandom[Finite] = Rand
		RandomIntersections.append(int(IntersectEctopicMatrices(EctopicArrayExp, EctopicPredRandom, sigma)))
	ExpPredIntersection = int(IntersectEctopicMatrices(EctopicArrayExp, EctopicArrayPred, sigma))
	return {
		"Random": json.dumps(RandomIntersections),
		"Real": ExpPredIntersection
		}

def EctopicGraphArray(RawArray):
	PlottingArray = numpy.copy(RawArray)
	PlottingArray = PlottingArray.astype('float')
	PlottingArray[numpy.logical_and(PlottingArray < 2, PlottingArray > -2)] = numpy.nan
	return json.dumps(PlottingArray[150:350, 150:350].tolist())

# ------======| DRAFT VISUALIZATION |======------

def VisualizeCool(InputCool, OutputPng, Region):
	subprocess.Popen(f"cooler show --out \"{OutputPng}\" -b --dpi 150 \"{InputCool}\" {Region}", shell = True, executable = "/bin/bash", stderr = subprocess.PIPE).communicate()

def VisualizePR(PRData, Name, FN):
	Precision, Recall = json.loads(PRData["Precision"]), json.loads(PRData["Recall"])
	fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
	ax.plot(Recall, Precision)        
	ax.set(xlabel="Recall", ylabel="Precision", title=f"{Name} PR Curve\nAUC = {PRData['AUC']:.10f}")
	ax.grid()
	fig.savefig(FN)
	plt.clf()

def VisualizeRandom(RandomData, FN):
	Random, Real = json.loads(RandomData["Random"]), int(RandomData["Real"])
	fig, ax = plt.subplots(figsize=(8,6), dpi=150)
	ax.set(title=f"Real vs Random Ectopic Intersections")
	ax.hist(Random, bins=200, histtype='step')
	ax.axvline(x=Real, color="red")
	fig.savefig(FN)
	plt.clf()

def VisualizeEctopicArray(EctopicArray, FN):
	Data = json.loads(EctopicArray)
	fig, ax = plt.subplots(figsize=(6, 6), dpi=150)
	ax.set(title=f"Ectopic Interactions")
	pos = ax.matshow(Data, cmap="bwr", vmin=-8, vmax=8)
	fig.colorbar(pos, ax=ax)
	fig.savefig(FN)
	plt.clf()

# ------======| HASH |======------

def HashJSON(Object, Key): return hashlib.sha256((json.dumps(Object, ensure_ascii=False) + Key).encode('utf8')).hexdigest()

def CheckHash(Object, key):
	Object = copy.deepcopy(Object)
	ObjectHash = copy.deepcopy(Object["__hash__"])
	del Object["__hash__"]
	return HashJSON(Object, key) == ObjectHash

# ------======| PARSER |======------

def CreateParser():
	Parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=f"BenchmarkPipeline by {__author__}", epilog=f"Email: regnveig@ya.ru")
	Parser.add_argument('--version', action='version', version=__version__)
	
	Parser.add_argument('-a', '--author', required=True, type=str, help=f"Author ID")
	Parser.add_argument('-l', '--model', required=True, type=str, help=f"Model name")
	Parser.add_argument('-s', '--sample', required=True, type=str, help=f"Sample name")
	Parser.add_argument('-r', '--resolution', required=True, type=int, help=f"Resolution")
	Parser.add_argument('-t', '--table', required=True, type=str, help=f"Main table path")
	Parser.add_argument('-w', '--wildtype', required=True, type=str, help=f"WT prediction TSV path")
	Parser.add_argument('-m', '--mutation', required=True, type=str, help=f"Mut prediction TSV path")
	Parser.add_argument('-o', '--output', required=True, type=str, help=f"Output dir")
	
	return Parser

# ------======| FILE CREATOR |======------

def CreateDataFiles(AuthorName, ModelName, SampleName, FileNamesInput, FileNamesOutput, Chrom, CaptureStart, CaptureEnd, RearrStart, RearrEnd, BinSize):
	
	# Create ts
	SubmissionDate = datetime.datetime.now(datetime.timezone.utc).isoformat()
	
	# Create data struct
	Data = { 
		"Metadata": {
			"AuthorID": AuthorName,
			"SubmissionDate": SubmissionDate,
			"ModelName": ModelName,
			"SampleName": SampleName,
			"Resolution": BinSize,
			"Chrom": Chrom,
			"Capture": { "Start": CaptureStart, "End": CaptureEnd },
			"Rearrangement": { "Start": RearrStart, "End": RearrEnd}
			}, 
		"Metrics": {}
		}
	
	# Create temp files
	with Timer(f"Temp files created") as _:
		TempDir = tempfile.TemporaryDirectory()
		TempFiles = {Type: os.path.join(TempDir.name, f"{'-'.join(Type)}.cool") for Type in FileNamesInput.keys()}
		for Type, FN in FileNamesInput.items():
			if Type[1] == "Exp": Cool2Cool(FN, TempFiles[Type], Chrom)
			if Type[1] == "Pred": Tsv2Cool(FN, TempFiles[Type], TempFiles[("Wt", "Exp")], Chrom, BinSize)
		Datasets = {Type: cooler.Cooler(FN) for Type, FN in TempFiles.items()}
	
	# Align
	with Timer(f"Sample type align") as _:
		SampleTypeAligned = {Type: os.path.join(TempDir.name, f"{'-'.join(Type)}-SampleTypeAligned.cool") for Type in FileNamesInput.keys()}
		for ST in ["Wt", "Mut"]: AlignCools(TempFiles[(ST, "Exp")], TempFiles[(ST, "Pred")], SampleTypeAligned[(ST, "Exp")], SampleTypeAligned[(ST, "Pred")], Chrom)
		SampleTypeAligned = {Type: cooler.Cooler(FN) for Type, FN in SampleTypeAligned.items()}
		
		# DRAFT
		for Type, Cool in SampleTypeAligned.items(): 
			VisualizeCool(InputCool = Cool.store, OutputPng = os.path.join(FileNamesOutput["DraftDir"], f"{Type[0]}{Type[1]}SampleTypeAligned.png"), Region = f"{Chrom}:{CaptureStart}-{CaptureEnd}")
		
	with Timer(f"Data type align") as _:
		DataTypeAligned = {Type: os.path.join(TempDir.name, f"{'-'.join(Type)}-DataTypeAligned.cool") for Type in FileNamesInput.keys()}
		for DT in ["Exp", "Pred"]: AlignCools(TempFiles[("Wt", DT)], TempFiles[("Mut", DT)], DataTypeAligned[("Wt", DT)], DataTypeAligned[("Mut", DT)], Chrom)
		DataTypeAligned = {Type: cooler.Cooler(FN) for Type, FN in DataTypeAligned.items()}
	
	# METRICS
	
	# Pearson
	with Timer(f"Pearson") as _:
		Data["Metrics"]["Pearson"] = {
			"WT": PearsonCorr(GetMatrix(SampleTypeAligned[("Wt", "Exp")], Chrom)["balanced"], GetMatrix(SampleTypeAligned[("Wt", "Pred")], Chrom)["balanced"]),
			"Mut": PearsonCorr(GetMatrix(SampleTypeAligned[("Mut", "Exp")], Chrom)["balanced"], GetMatrix(SampleTypeAligned[("Mut", "Pred")], Chrom)["balanced"]),
			}

	# SCC
	with Timer(f"SCC") as _:
		Data["Metrics"]["SCC"] = {
			"WT": SCC(SampleTypeAligned[("Wt", "Exp")], SampleTypeAligned[("Wt", "Pred")], MaxDist = C_SCC_MAXDIST, h = C_SCC_H),
			"Mut": SCC(SampleTypeAligned[("Mut", "Exp")], SampleTypeAligned[("Mut", "Pred")], MaxDist = C_SCC_MAXDIST, h = C_SCC_H)
			}

	# Insulation
	with Timer(f"Insulation Dataset") as _:
		InsDataset = InsulationData(Datasets, Window=BinSize * 5)
	
	with Timer(f"Insulation Score Pearson") as _:
		Data["Metrics"]["InsulationScorePearson"] = {
			"WT": PearsonCorr(InsDataset["sum_balanced_Wt-Exp"], InsDataset["sum_balanced_Wt-Pred"]),
			"Mut": PearsonCorr(InsDataset["sum_balanced_Mut-Exp"], InsDataset["sum_balanced_Mut-Pred"])
			}
	
	with Timer(f"Insulation Score (Mut/Wt) Pearson") as _:
		Data["Metrics"]["InsulationScoreMutVsWtPearson"] = PearsonCorr(InsDataset["sum_balanced_Mut/Wt-Exp"], InsDataset["sum_balanced_Mut/Wt-Pred"])

	# Insulatory AUC
	with Timer(f"Ectopic Insulation PR Curve") as _:
		Data["Metrics"]["EctopicInsulationPR"] = PRCurve(
			YTrue = InsDataset["Y-True"],
			Probas = InsDataset["sum_balanced_sigma_Mut/Wt-Pred"]
			)
		
		# DRAFT
		VisualizePR(PRData = Data["Metrics"]["EctopicInsulationPR"], Name = "Ectopic Insulation", FN = os.path.join(FileNamesOutput["DraftDir"], f"EctopicInsulationPRCurve.png"))

	# Ectopic
	with Timer(f"Ectopic Array") as _:
		EctopicArrayExp = EctopicInteractionsArray(Datasets[("Wt", "Exp")], Datasets[("Mut", "Exp")], Chrom, CaptureStart, CaptureEnd, RearrStart, RearrEnd, Normalized=False)
		EctopicArrayPred = EctopicInteractionsArray(Datasets[("Wt", "Pred")], Datasets[("Mut", "Pred")], Chrom, CaptureStart, CaptureEnd, RearrStart, RearrEnd, Normalized=True)
	
	with Timer(f"Ectopic Array Graph") as _:
		Data["Metrics"]["EctopicArrayGraph"] = {
			"Exp": EctopicGraphArray(EctopicArrayExp),
			"Pred": EctopicGraphArray(EctopicArrayPred)
			}
		
		# DRAFT
		for Key, EctopicArray in Data["Metrics"]["EctopicArrayGraph"].items(): VisualizeEctopicArray(EctopicArray = EctopicArray, FN = os.path.join(FileNamesOutput["DraftDir"], f"{Key}EctopicArray.png"))
	
	with Timer(f"Ectopic Interactions PR Curve") as _:
		Data["Metrics"]["EctopicInteractions"] = PRCurve(
			YTrue = [(i == i) and ((i > 2) or (i < -2)) for i in EctopicArrayExp.flatten()],
			Probas = EctopicArrayPred.flatten()
			)
		
		# DRAFT
		VisualizePR(PRData = Data["Metrics"]["EctopicInteractions"], Name = "Ectopic Interactions", FN = os.path.join(FileNamesOutput["DraftDir"], f"EctopicInteractionsPRCurve.png"))
	
	# Random
	with Timer(f"Random Interactions") as _:
		Data["Metrics"]["RandomInteractions"] = RandomEctopicIntersections(EctopicArrayExp, EctopicArrayPred)
		
		# DRAFT
		VisualizeRandom(RandomData = Data["Metrics"]["RandomInteractions"], FN = os.path.join(FileNamesOutput["DraftDir"], f"RealVsRandomEctopicInteractions.png"))
	
	# HASH
	with Timer(f"HASH") as _:
		Data["__hash__"] = HashJSON(Data, "lol")
	
	# SAVE
	with Timer(f"Save JSON") as _:
		json.dump(Data, open(FileNamesOutput["Json"], 'wt'), ensure_ascii=False, indent=4)
	
	with Timer(f"Save MCOOLs") as _:
		for Key in SampleTypeAligned.keys(): MakeMcool(InputCool = SampleTypeAligned[Key].store, OutputMcool = FileNamesOutput[Key], Resolution = BinSize)

# ------======| MAIN |======------

def Main():
	
	Parser = CreateParser()
	Namespace = Parser.parse_args(sys.argv[1:])
	
	ConfigureLogger("benchmark_pipeline.log")
	
	# Load Main Table
	MainTable = pandas.read_csv(Namespace.table, sep='\t', dtype=str).set_index("rearrangement_ID").rename_axis(None, axis=0).transpose().to_dict()
	try:
		SampleData = MainTable[Namespace.sample]
	except KeyError:
		raise ValueError(f"Unknown Sample Name: '{Namespace.sample}'")
	
	FileNamesInput = {
		("Wt", "Exp"): os.path.join(SampleData["capture_WT_raw_data"], f"inter_{int(Namespace.resolution / 1000)}kb.cool"),
		("Wt", "Pred"): Namespace.wildtype,
		("Mut", "Exp"): os.path.join(SampleData["capture_Mut_raw_data"], f"inter_{int(Namespace.resolution / 1000)}kb.cool"),
		("Mut", "Pred"): Namespace.mutation
		}

	FileNamesOutput = {
		"Json": os.path.join(Namespace.output, "metadata.json"),
		"DraftDir": os.path.join(Namespace.output, ".draft"),
		("Wt", "Exp"): os.path.join(Namespace.output, "WtExp.cool"),
		("Wt", "Pred"): os.path.join(Namespace.output, "WtPred.cool"),
		("Mut", "Exp"): os.path.join(Namespace.output, "MutExp.cool"),
		("Mut", "Pred"): os.path.join(Namespace.output, "MutPred.cool")
		}
	
	os.mkdir(Namespace.output)
	os.mkdir(FileNamesOutput["DraftDir"])
	
	CreateDataFiles(
		AuthorName = Namespace.author,
		ModelName = Namespace.model,
		SampleName = Namespace.sample,
		FileNamesInput = FileNamesInput, 
		FileNamesOutput = FileNamesOutput,
		Chrom = SampleData["chr"],
		CaptureStart = int(SampleData["start_capture"]),
		CaptureEnd = int(SampleData["end_capture"]),
		RearrStart = int(SampleData["start1"]),
		RearrEnd = int(SampleData["end1"]),
		BinSize = Namespace.resolution)

if __name__ == "__main__": Main()
