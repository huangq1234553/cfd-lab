#! /usr/bin/env python3

### Rescales .pgm geometry files
# Note: forbidden geometries will NOT be checked/corrected! Also there is no guarantee
#   that rescaling will not introduce any forbidden geometry!
#
# Before usage you might need to make this file executable with:
#   chmod +x rescale-geometry.sh
#
# Usage:
#   ./rescale-geometry.py SCALE_FACTOR INPUT_FILE OUTPUT_FILE
#

from sys import argv, exit

###
class InOutFileTranslator:
    _inFile = None
    _outFile = None
    _linesQ = []
    _safeLinesCount = 0
    _translateFunction = None

    def __init__(self, inputFileName, outputFileName, translateFunction):
        """
        Initialize the translator
        :param inputFileName: the path to input file
        :param outputFileName: the path to output file
        :param translateFunction: a function accepting a single string argument which can be used to translate lines
        """
        self._inFile = open(inputFileName, 'r')
        self._outFile = open(outputFileName, 'w')
        self._initializeLineQ(self._inFile.readlines())
        self._inFile.close()
        self._translateFunction = translateFunction

    def close(self):
        self._outFile.close()

    def isRawLineSafe(self, line):
        stripped = line.strip()
        return not( len(stripped) == 0 or stripped[0] == "#" )

    def isLineSafe(self, lineTuple):
        return lineTuple[0]

    def _initializeLineQ(self, rawList):
        for line in rawList:
            if self.isRawLineSafe(line):
                self._linesQ.append((True, line))
                self._safeLinesCount += 1
            else:
                self._linesQ.append((False, line))
        self._linesQ.reverse() # so we can pop from file's head

    def rawCount(self):
        return len(self._linesQ)

    def count(self):
        return self._safeLinesCount

    def read(self):
        """
        Here just read next line and return it if not a comment (or empty line). If it is a comment, just write it, as is, to the output file.
        :return: string containing the next non-comment line, IndexError is thrown in case file's end is reached
        """
        cur = self._linesQ.pop()
        while not self.isLineSafe(cur):
            self._outFile.write(cur[1])
            cur = self._linesQ.pop()
        self._safeLinesCount -= 1
        return cur[1]

    def write(self, line):
        self._outFile.write(line)

    def translateLine(self, count=1):
        """
        Translates current input line and writes it to output, additionally it returns a copy of the translated string.
        :param count: the number of copies of the translated line to be written to output
        :return: A copy of translated string.
        """
        cur = self._translateFunction(self.read())
        for i in range(count):
            self.write(cur)
        return cur

def managePgmHeader(translator, factor):
    # Here just parse the header and make sure to adapt it into output file
    magic = translator.read()
    translator.write(magic)
    dimensions = translator.read()
    W,H = [ int(x) for x in dimensions.strip().split(' ') ]
    nW = (W-2)*factor + 2
    nH = (H-2)*factor + 2
    translator.write("%d %d\n" %(nW,nH))
    depth = translator.read()
    translator.write(depth)

def rescaleLine(line, factor):
    lineElements = line.strip().split(' ') # be careful, strip also removes the newline char at the end!
    oLine = [lineElements[0]]
    for x in lineElements[1:-1]:
        for i in range(factor):
            oLine.append(x)
    oLine.append(lineElements[-1]) # Here we just append the boundary ...
    oLine.append('\n') # ... and the newline
    return ' '.join(oLine)

def rescaleImage(translator, factor):
    translator.translateLine()  # Top boundary
    while translator.count() > 1:
        translator.translateLine(count=factor)
    translator.translateLine() # Bot boundary

def maskLine(line, valuesTo1):
    lineElements = line.strip().split(' ') # be careful, strip also removes the newline char at the end!
    oLine = []
    for x in lineElements:
        if x in valuesTo1:
            x = '1'
        else:
            x = '0'
        oLine.append(x)
    oLine.append('\n') # ... and the newline
    return ' '.join(oLine)
###

USAGE="Usage:\n\t./rescale-geometry.py SCALE_FACTOR INPUT_FILE OUTPUT_FILE [OUTPUT_MASK_FILE]\nAll files must be .pgm!"

if (len(argv) == 4):
    scaleFactor, inputFileName, outputFileName = argv[1:]
    outputMaskFileName = None
elif (len(argv) == 5):
    scaleFactor, inputFileName, outputFileName, outputMaskFileName = argv[1:]
else:
    print(USAGE)
    exit(1)

scaleFactor = int(scaleFactor)

translateFunction = lambda x: rescaleLine(x, scaleFactor)

pgmTranslator = InOutFileTranslator(inputFileName, outputFileName, translateFunction)

managePgmHeader(pgmTranslator, scaleFactor)
rescaleImage(pgmTranslator, scaleFactor)

pgmTranslator.close()

if outputMaskFileName:
    maskFunction = lambda x: maskLine(x, ['0', '1']) # we want to mask all no-slip and free-slip cells by default
    maskTranslator = InOutFileTranslator(outputFileName, outputMaskFileName, maskFunction)
    managePgmHeader(maskTranslator, 1)
    rescaleImage(maskTranslator, 1)
    maskTranslator.close()

#eof
