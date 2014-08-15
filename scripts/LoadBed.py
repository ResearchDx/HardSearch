import os
import collections

def LoadBedFile(BedFilePath):
    realBedPath = os.path.realpath(BedFilePath)
    resultDict = collections.defaultdict(list)
    try:
        bedFile = open(realBedPath, 'r')
        line = bedFile.readline()
        lineCount = 0
        while line:
            data = line.split('\t')
            if data[0] != 'Chr':
                resultDict[data[0].lower()].append((data[1], data[2]))
            line = bedFile.readline()
    except IOError:
        print("Bed file not found, target regions empty")
        return None
    else:
        return resultDict

def IsIntersecting(pos, bedCoordinates):
    ''' Returns true if position is between the two coordinate tuple, otherwise
        false
    '''
    if pos >= int(bedCoordinates[0]) and pos <= int(bedCoordinates[1]):
        return True
    return False

def IsIntersectingBed(pos, bedDict):
    for bedCoordinate in bedDict[pos[0]]:
        if IsIntersecting(pos[1], bedCoordinate):
                return True
        elif IsIntersecting(pos[2], bedCoordinate):
                return True
        elif IsIntersecting(int(bedCoordinate[0]), (pos[1], pos[2])):
                return True
        elif IsIntersecting(int(bedCoordinate[1]), (pos[1], pos[2])):
                return True
    return False
