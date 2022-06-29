""" ******************************************

            Assignment 1
Name: Jnaneswara Rao Rompilli (EE20B052)
Date: 24-01-2022
Description: Take circuit parameters as input and extract circuit parameters
Output: Reversed values and lines

****************************************** """

from sys import argv, exit


CIRCUIT = ".circuit"
END = ".end"


# Extract the required data from lines (remove comments)
def extractData(lines_raw):
    lines_actual = []
    for i in range(0, len(lines_raw)):
        line = lines_raw[i].split("#")[0].split()  # remove comments and trailing spaces
        lines_actual.append(" ".join(line))

    return lines_actual


# Extract tokens from a Line
def extractTokens(Line):
    words = Line.split()

    if len(words) == 4:
        elementName, node1, node2, value = [x for x in words]
        return [elementName, node1, node2, value]

    # CCxS
    elif len(words) == 5:
        elementName, node1, node2, voltageSource, value = [x for x in words]
        return [elementName, node1, node2, voltageSource, value]

    # VCxS
    elif len(words) == 6:
        elementName, node1, node2, voltageSourceN1, voltageSourceN2, value = [x for x in words]
        return [elementName, node1, node2, voltageSourceN1, voltageSourceN2, value]

    else:
        return []


# Print values in mentioned format
def printCircuit(lines):
    output = ""
    for i in reversed(range(0, len(lines))):  # iterate over valid range
        line = lines[i].split()  # split line into words
        line.reverse()  # reverse the list
        output = output + (" ".join(line) + "\n")  # join words after reversing and add "\n" at the end of line

    print(output, end="\0")


if len(argv) != 2:  # If more than 2 arguments are passed
    print("Invalid operation !")
    print(f"Usage: {argv[0]} <inputfile>'")
    exit()

try:
    with open(argv[1]) as f:
        lines = f.readlines()
        start = -1
        end = -2
        start_cnt = 0
        end_cnt = 0
        for line in lines:  # extracting circuit definition: start and end lines
            if CIRCUIT == line[0 : len(CIRCUIT)]:
                start = lines.index(line)
                start_cnt = start_cnt + 1
            elif END == line[: len(END)]:
                end = lines.index(line)
                end_cnt = end_cnt + 1

        if (start_cnt > 1 or end_cnt > 1):  # Check if there are multiple .circuit/.end declarations in input file
            print("Invalid circuit definition!\nMultiple .circuit/.end declarations detected")
            exit()

        if start >= end:  # validating circuit block
            print("Invalid circuit definition!")
            exit()

        LinesRaw = lines[start + 1 : end]
        LinesActual = extractData(LinesRaw)
        LinesTokens = [extractTokens(line) for line in LinesActual]

        printCircuit(LinesActual)
        # print(LinesTokens)


except IOError:
    print("Invalid file!")
    exit()
