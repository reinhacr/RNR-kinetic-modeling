import os
import re
import shutil
from glob import glob
import subprocess
from time import time, sleep
import numpy as np
import csv
import sys
import argparse
from contextlib import redirect_stderr
import pickle
from itertools import zip_longest


def readSimpleInput(paramFile, getArgs=None):
    """read input file with the simple form:
    param1: value
    param2: value, etc.

    Arguments
    ---------
    paramFile: str
        any file containing input parameters in correct form
    args: argparse.Namespace
        existing parsed args to modify (optional)

    Returns
    -------
    args: argparse.Namespace
        either modified original or new
    """

    try:
        args = argparse.Namespace()
        with open(paramFile, "r") as infile:
            text = infile.read()
            text = text.split("\n")
            text = list(filter(None, text))

            for row in text:
                row = row.replace("\t", "")
                row = row.strip("\n")
                row = row.lstrip(" ")
                row = row.rstrip(" ")
                if row[0] == "#":  # ignore comment lines
                    continue
                key = row.split(":")[0].strip(" ")
                value = row.split(":")[1].strip(" ")
                if value.lower() == "true":
                    value = True
                elif value.lower() == "false":
                    value = False
                elif re.search(r"^-?([0-9])+$", value):
                    value = int(value)
                elif re.search(r"^-?([0-9]*)(\.)([0-9])*$", value):
                    value = float(value)
                elif value.lower() == "none":
                    value = None
                elif "," in value:
                    # list
                    value = value.replace(" ", "")
                    value = value.split(",")
                    clean = []
                    for item in value:
                        if item.lower() == "true":
                            clean.append(True)
                        elif item.lower() == "false":
                            clean.append(False)
                        elif re.search(r"^[0-9][\.]?[0-9]*$", item):
                            clean.append(float(item))
                        elif item.lower() == "none":
                            clean.append(None)
                        else:
                            clean.append(item)
                    value = clean

                setattr(args, key, value)

        if getArgs is not None:
            try:
                if type(getArgs) == str:
                    return vars(args)[getArgs]
                else:
                    argList = []
                    for arg in getArgs:
                        argList.append(vars(args)[arg])
                    return argList
            except KeyError:
                print("input file missing requested parameter(s)")
        else:
            return args

    except FileNotFoundError:
        sys.exit("missing input file")

    return args


def cyCompile(script, opt_level=3, name=None):
    """automates simple Cythonizing, called by cyCompile command line tool"""

    name = "".join(script.split(".pyx")[0:-1]) if name is None else name

    with open(f"{name}_setup.py", "w+") as setup:
        setup.write(f"""
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from os import environ, system as sh
from shutil import rmtree
import numpy as np

# environ["CC"] = "/usr/bin/gcc"
extension = Extension(
    name="{name}",
    sources=["{script}"],
    include_dirs=["/usr/bin/include", np.get_include()],
    library_dirs=["/usr/bin/lib"],
    libraries=['m'],
    extra_compile_args=['-O3', '-march=x86-64'],
    extra_link_args=['-O3', '-march=x86-64'])
setup(name="{name}", ext_modules=cythonize([extension], build_dir="."))
sh("mv {name}.cpython* {name}.so")
sh("rm {name}.c")
rmtree("./build")""")

    out = os.system(f"python {name}_setup.py build_ext --inplace &>{name}.log")
    if out != 0:
        print(f"Cython build failed, see {name}.log")
    else:
        print(f"{name} built successfully")


# class flexParse
def flexParse(ap: argparse.ArgumentParser):
    """untested!!!"""
    try:
        with open(os.devnull, "w") as void:
            with redirect_stderr(void):
                args = ap.parse_args()
    except SystemExit:
        try:
            newAp = argparse.ArgumentParser()
            newAp.add_argument("infile")
            with open(os.devnull, "w") as void:
                with redirect_stderr(void):
                    infileArg = newAp.parse_args()

            moreArgs = readInputFile(infileArg.infile)
            argsDict = vars(args[0])  # remove unfilled args list
            argsDict.update(vars(moreArgs))
            args = argparse.Namespace(**argsDict)

        except SystemExit or FileNotFoundError:
            print("need an infile or arguments here")

    return args


def combineResults(outfile, average=False, ext="results_temp"):
    """average or concatenate data from parallel processes
    data must be in csv form and in columns"""

    if average is True:
        data = []
        for datafile in glob(f"*.{ext}"):
            nHeaders = 0
            with open(datafile, "r") as result:
                line = result.readline()
                # scroll until hit numbers
                # find last non-numeric line for column names, if present
                while (re.match("^[A-Za-z].*$", line)) and line:
                    line = result.readline()
                    nHeaders += 1
            currData = np.genfromtxt(
                datafile, delimiter=",", skip_header=nHeaders)
            # nan occurs if headers have more cols than data
            currData = [row[~np.isnan(row)] for row in currData]
            data.append(currData)

        data = np.array(data)
        if type(data[0]) == list:
            # happens when subdata has different sizes
            raise ValueError("different lengths")
        else:
            avgData = np.mean(data, axis=0)

        with open(outfile, "w+") as avgResults:
            writer = csv.writer(avgResults)
            for line in avgData:
                writer.writerow(line)
    else:
        # concatenate in correct order
        files = glob(f"*.{ext}")
        idx = []
        for datafile in files:
            idx.append(int(datafile.split(".")[0]))

        idx = sorted(idx)
        files = [f"{idx[i]}.{ext}" for i in range(len(files))]

        # binary I/O is faster
        with open(outfile, "wb+") as out:
            for datafile in files:
                print(datafile)
                with open(datafile, "rb") as curr:
                    data = curr.read()
                    # linebreak is extremely important!
                    if data[-1] != b"\n":
                        data = data + b"\n"
                    out.write(data)


def readInputFile(fileName, getArg=None):
    """parse files like example.input"""

    args = DataTree()

    with open(fileName, "r") as infile:
        lines = infile.readlines()
        lines = [lines[i] for i in range(len(lines))
                 if not lines[i].isspace()]
        row = 0
        # find beginning of parameters
        while not lines[row].startswith("ARGUMENTS"):
            row += 1
        else:
            row += 2

        # get arguments
        lines = lines[row:]
        lines = [line for line in lines if "_|" in line]

        for row in range(len(lines)):
            key = lines[row].split("_")[0]
            lines[row] = lines[row].split("| ")[1]
            lines[row] = lines[row].split("#")[0]  # comments
            arg = lines[row].strip()  # residual spaces, tabs, newlines
            if re.search("^[0-9]+$", arg):
                arg = int(arg)
            elif re.search("^[0-9]+\.[0-9]+$", arg):
                arg = float(arg)
            elif re.search("[0-9]*[\.]?[0-9]*[eE][0-9]+", arg):
                arg = float(arg)
            elif (arg.startswith("[") and arg.endswith("]")) or\
                 (arg.startswith("(") and arg.endswith(")")):
                arg = arg[1:-1]
                arg = re.sub("[\, ]", "|", arg)
                arg = list(filter(None, arg.split("|")))
                # now try to convert to numbers
                for i in range(len(arg)):
                    try:
                        arg[i] = float(arg[i])
                    except ValueError:
                        pass
            elif arg.lower() in ["true", "false"]:
                arg = True if arg.lower() == "true" else False
            elif arg.lower() == "none":
                arg = None
            args.add(key, arg)

        # adds numbers to outfile to avoid overwrite

        args.out = overwriteProtect(args.out)

    if getArg is not None:
        try:
            if type(getArg) == str:
                return vars(args)[getArg]
            else:
                argList = []
                for arg in getArg:
                    argList.append(vars(args)[arg])
                return argList
        except KeyError:
            print("input file missing requested parameter(s)")
    else:
        return args


def padArray(arr, withWhat=np.nan):
    """even out list of lists or irregular 2D ndarray with for
    later use with NumPy's fast methods or easy printing"""

    if (type(arr) == list) and (len(arr[0]) < 2):
        raise ValueError("array must be 2D")
    elif (type(arr) == np.ndarray) and (len(arr.shape) != 2):
        raise ValueError("array must be 2D")
    else:
        lens = [len(row) for row in arr]
        maxLen = max(lens)

        padded = []
        for row, length in zip(arr, lens):
            padded.append(list(row) + [withWhat] * (maxLen - length))

        return np.array(padded)


def findMissing(fileList: str, direc: str="./") -> None:
    if not direc.endswith("/"):
        direc = direc + "/"
    missing = []
    for file in fileList:
        if not os.path.exists(direc + file):
            missing.append(file)
    if missing:
        raise RuntimeError(f"missing source files: {missing}")


def waitFor(func, args=None, corr=True, checkEvery=1):
    """emulates "until" in bash scripts"""
    if args is None:
        val = func()
        while val != corr:
            sleep(checkEvery)
            val = func()
    else:
        val = func(*args)
        while val != corr:
            sleep(checkEvery)
            val = func(*args)


def rm(string: str) -> None:
    patterns = string.split(" ")
    for pattern in patterns:
        for file in glob(pattern):
            os.remove(file)


def rmDir(string: str) -> None:
    patterns = string.split(" ")
    for pattern in patterns:
        for direc in glob(pattern):
            shutil.rmtree(direc)


def mv(string: str, dest: str) -> None:
    strings = string.split(" ")
    for string in strings:
        for file in glob(string):
            shutil.move(file, dest)


def slurm(script: str) -> int:
    try:
        out = os.system(f"sbatch {script} &>/dev/null")
    except FileNotFoundError:
        out = os.system(f"sbatch {glob(script)[0]} &>/dev/null")
    if out != 0:
        sys.exit("slurm submission failed")

    return out


def submitJobArray(arrayFile: str, partition: str, cores: int,
                   time: int, memGB: int=5) -> None:

    memGB = int(memGB)  # just in case
    os.system("conda deactivate &>/dev/null; module load dSQ; " +
              f"dsq --job-file {arrayFile} -p {partition} -N 1 -n 1" +
              f" -c {cores} -t {time}:00 --mem-per-cpu={memGB}GB " +
              "--requeue &>/dev/null;")
    out = slurm("dsq*.sh")
    if out != 0:
        sys.exit("dsq submission failed, file does not exist")


def slurmCount(jobname, status=None):
    """
    Arguments
    ---------
    jobname: str
    status: str
        "running", "pending", defaults to "all"

    Returns
    -------
        nJobs: int
    """

    # checks the length of the *.tsv job file instead of calling squeue
    # to avoid slowing down Grace for other people
    # nJobs = subprocess.run('grep -c ".*" *.tsv', shell=True, stdout=subprocess.PIPE).stdout

    if status is None:
        nJobs = subprocess.run(
            "read -ra arr <<< " +
            f"$(squeue --name {jobname} -u $(whoami)); " +
            "len=${#arr[@]}; n=$(($((len / 11)) - 1)); echo $n;",
            shell=True, stdout=subprocess.PIPE).stdout
    else:
        nJobs = subprocess.run(
            "read -ra arr <<< " +
            f"$(squeue --name {jobname} -u $(whoami) -t {status}); " +
            "len=${#arr[@]}; n=$(($((len / 11)) - 1)); echo $n;",
            shell=True, stdout=subprocess.PIPE).stdout

    return int(nJobs.strip())


def mkdir(string: str) -> None:
    for direc in string.split(" "):
        try:
            os.mkdir(direc)
        except FileExistsError:
            pass


def addEnv(**kwargs):
    """emulates 'export' cmd in shell"""
    for key, value in kwargs.items():
        try:
            original = os.environ[key]
            os.environ[key] = str(value) + ":" + original
        except KeyError:
            os.environ[key] = str(value)


def sendEmail(*, subject, content, recipient):
    """wrapper for shell mail"""

    if os.path.exists(content):
        os.system(f'mail -s "{subject}" {recipient} < {content}')
    else:
        # check this
        os.system(f'mail -s "{subject}" {recipient} <<< {content}')


def exists(filepath):
    if os.path.exists(filepath):
        return True
    else:
        return False


def direcs(pattern="*", path="."):
    """ equivalent of $(ls -p | grep / | grep [pattern]) """
    direcs = []
    path = path + "/" if path[-1] != "/" else path
    for thing in os.listdir(path):
        fullPath = path + thing
        if os.path.isdir(fullPath) and thing in glob(pattern):
            direcs.append(thing)
    return direcs


def replaceText(filename, pattern, replacement, newFile=None):
    """mimics sed, but slightly different

    EXAMPLE: replaceText("data.txt", "elephant", 'frog")
    this replaces every instance of 'elephant' with 'frog'

    EXAMPLE: replaceText("data.txt", r"man.*\\b::man", "frog")
    this finds every word that STARTS with 'man' and replaces
    the string 'man' within these words with 'frog'. The "::"
    separates the search pattern from the part to be replaced
    """

    with open(filename, "r") as f:
        text = f.read()

    if "::" in pattern:
        try:
            searchPattern, pattern = pattern.split("::")
        except ValueError:
            sys.exit("cannot use more than one '::' in pattern")

        hits = re.findall(searchPattern, text)
        replaced = [re.sub(pattern, replacement, hit) for hit in hits]
        for hit in hits:
            text = re.sub(hit, replaced, text)
    else:
        text = re.sub(pattern, replacement, text)

    if newFile is not None:
        with open(newFile, "w+") as f:
            f.write(text)
    else:
        with open(filename, "w+") as f:
            f.write(text)


def overwriteProtect(filename):
    """adds _1, _2, etc. to filename if name is taken to avoid overwriting"""

    if filename in os.listdir():
        stem = filename.split(".")[0]
        ext = filename.split(".")[-1]
        names = re.findall(f"{stem}_[0-9]+\.{ext}", ";;".join(os.listdir()))
        # find highest number
        highestNum = 0
        for name in names:
            newNum = int(re.findall(f"[0-9]+", name)[0])
            if newNum > highestNum:
                highestNum = newNum

        return f"{stem}_{highestNum + 1}.{ext}"
    else:
        return filename


class DataTree:
    """
    Make, read, and write arbitrary data structures with ease

    Methods
    -------
    add(key, value=None)
        add a new leaf to the current branch
    remove(key)
        remove a leaf (not a branch)
    newBranch(name)
        create a new DataTree object under the current one
    removeBranch(name)
        remove an entire branch
    load(filename)
        load a pickled DataTree from file to the current Tree
    write(filename)
        save DataTree in reconstructable binary form
    load_csv(filename)
        load column-wise numerical data with column headers from csv
    write_csv(filename)
        write column-wise data with column headers to csv
    """

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def __str__(self):
        heading = "DataTree: \n---------"

        def unpack(data, level=0, structure=[]):
            for key in data:
                # indent according to level in heirarchy
                structure.append(level * "    " + "- " + key)

                subStruct = data[key]
                if type(subStruct) == dict:
                    # move down a level recursively
                    level += 1
                    unpack(data[key], level, structure)
                    # return to previous level
                    level -= 1
                else:
                    try:
                        # find length if any and print it
                        length = subStruct.__len__()
                        structure[-1] = structure[-1] + ": " +\
                            type(subStruct).__name__ + " (" + str(length) + ")"
                    except AttributeError:
                        # just print type without length
                        structure[-1] = structure[-1] + ": " + type(subStruct).__name__

            return structure

        structure = "\n".join(unpack(self.__dict__))
        return heading + "\n" + structure

    def __repr__(self):
        return self.__str__()

    def add(self, key, value=None):
        if isinstance(value, DataTree):
            raise TypeError("to add a branch, use DataTree.newBranch()")
        else:
            self.__dict__[key] = value

    def newBranch(self, name):
        self.__dict__[name] = DataTree()
        return self.__dict__[name]

    def remove(self, name):
        if isinstance(self.__dict__.name, DataTree):
            raise TypeError("to remove a branch, use DataTree.removeBranch()")
        else:
            del self.__dict__[name]

    def removeBranch(self, name):
        if not isinstance(self.__dict__[name], DataTree):
            raise TypeError(f"Attribute '{name}' is not a branch")
        else:
            del self.__dict__[name]

    def load(self, filename):
        try:
            self.__dict__ = pickle.load(filename).__dict__
        except FileNotFoundError:
            print(f"file {filename} not found")

    def write(self, filename):
        pickle.dump(self, filename)

    def load_csv(self, filename):
        with open(filename, "r") as infile:
            reader = csv.reader(infile)
            line = next(reader)
            while (line == []) and line:
                line = next(reader)
            headers = list(filter(None, line))

            try:
                float(headers[0])
                raise ValueError("data file missing column header(s)")
            except ValueError:
                pass

            headerSet = set()
            for elem in headers:
                if elem in headerSet:
                    raise ValueError("header column names have duplicates")
                else:
                    headerSet.add(elem)

            data = []
            line = next(reader)
            while True:
                try:
                    data.append(line)
                    line = next(reader)
                except StopIteration:
                    break

        data = list(zip(*data))
        tree = {}
        try:
            for i in range(len(data)):
                data[i] = list(filter(None, data[i]))
                data[i] = list(map(float, data[i]))
                tree[headers[i]] = np.array(data[i])
        except ValueError:
            print("Datatree.load_csv(): columns cannot contain non-numeric data")
        except IndexError:
                    print("Datatree.load_csv(): number of headers < number of columns")

    def write_csv(self, filename, **kwargs):
        with open(filename, "w+") as out:
            writer = csv.writer(out)
            writer.writerow(list(kwargs.keys()))
            data = zip_longest(*list(kwargs.values()))
            for row in data:
                writer.writerow(row)

    def show_structure(self):
        pass


def timeMe(func):
    """convenient timer decorator"""
    def wrapper(*args, **kwargs):
        start = time()
        result = func(*args, **kwargs)
        timer = time() - start
        print(f"\n>>> Runtime, {func.__name__}: {timer:.3g} sec")
        return result
    return wrapper


class TimeSmart:
    """various timers with user-friendly output

    timer = TimeSmart()
    timer.tic()         # start
    timer.toc()         # display runtime in comfortable units
    looptime()          # estimate remaining time in for loop

    example:
    timer = TimeSmart()
    for i in range(nTotal):
        ...[your_code]...
        timer.loopTime(nTotal, i)
    """

    import time as time

    def __init__(self):
        self.tPerIter = 0

    def tic(self):
        self.now = self.time.time()

    def toc(self):
        t = self.time.time() - self.now

        if t < 1e-3:
            T = t * 1e6
            unit = "microseconds"
        elif 1e-3 <= t < 1:
            T = t * 1000
            unit = "milliseconds"
        elif 1 <= t < 60:
            T = t
            unit = "seconds"
        elif 60 <= t < 3600:
            T = t / 60.0
            unit = "minutes"
        elif 3600 <= t < 86400:
            T = t / 3600
            unit = "hours"
        else:
            T = t / 86400
            unit = "days"

        return f"Runtime: {T:0.2f} {unit}"

    def loopTime(self, nTotal: int, nDone: int) -> str:
        t = self.time.time() - self.now

        tCurr = (1 / nDone) * ((nDone - 1) * self.tPerIter + t)
        tLeft = (nTotal - nDone) * tCurr
        self.tPerIter = tCurr

        if tLeft < 10:
            return "finishing up..."
        elif 10 <= tLeft < 60:
            T = tLeft
            unit = "sec"
        elif 60 <= tLeft < 3600:
            T = tLeft / 60.0
            unit = "min"
        elif 3600 <= tLeft < 86400:
            T = tLeft / 3600
            unit = "hrs"
        else:
            T = tLeft / 86400
            unit = "days"

        return f"{T:0.2f} {unit} left"
