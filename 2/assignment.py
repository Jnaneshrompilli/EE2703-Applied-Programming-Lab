""" ******************************************

            Assignment 2
Name: Jnaneswara Rao Rompilli (EE20B052)
Date: 24-01-2022

Description: This code can solve AC as well as DC circuits
             Supported components:
             1) Resistor
             2) Capacitor
             3) Inductor
             4) Voltage Source
             5) Current Source
             6) Voltage controlled Voltage source
             7) Voltage controlled Current source

Output: All node voltages and current through voltage sources

****************************************** """


from sys import argv, exit
import numpy as np
from math import sin, cos, radians

# Intializing Variables
CIRCUIT = ".circuit"
END = ".end"
AC = ".ac"

PI = 3.1415

R = "Resistor"
C = "Capacitor"
L = "Inductor"
V = "VoltageSource"
I = "CurrentSource"
E = "VCVS"
G = "VCCS"
F = "CCCS"
H = "CCVS"

ac_simulation = False

G_matrix = np.zeros((0, 0), dtype="complex")  # Conductance matrix
V_vector = np.zeros((1,), dtype="complex")  # Variable vector
I_vector = np.zeros((1,), dtype="complex")  # Vector of independent sources

rows = 0
columns = 0

frequency = float(0)
omega = float(0)

node_count = 0
VI_count = [0, 0]
VS_count, IC_count = 0, 0

# Class for resistor
class Resistor:
    # Constructor to intialize the object with input values
    def __init__(self, name, node1, node2, value):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.value = value

    # Function to modify the conducatnce matrix according to resistor in Modified nodal analysis
    def fillMatrix(self, G_matrix, nodes):
        k = nodes[self.node1]
        l = nodes[self.node2]

        if self.node1 == "GND":
            G_matrix[l][l] += 1.0 / self.value

        elif self.node2 == "GND":
            G_matrix[k][k] += 1.0 / self.value

        else:
            G_matrix[k][k] += 1.0 / self.value
            G_matrix[k][l] += -1.0 / self.value
            G_matrix[l][k] += -1.0 / self.value
            G_matrix[l][l] += 1.0 / self.value


# Class for capacitor
class Capacitor:
    # Constructor to intialize the object with input values
    def __init__(self, name, node1, node2, value):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.value = value

    # Function to modify the conducatnce matrix according to capacitor in Modified nodal analysis
    def fillMatrix(self, G_matrix, nodes):
        k = nodes[self.node1]
        l = nodes[self.node2]

        if self.node1 == "GND":
            G_matrix[l][l] += 1.0 / self.value

        elif self.node2 == "GND":
            G_matrix[k][k] += 1.0 / self.value

        else:
            G_matrix[k][k] += 1.0 / self.value
            G_matrix[k][l] += -1.0 / self.value
            G_matrix[l][k] += -1.0 / self.value
            G_matrix[l][l] += 1.0 / self.value


# Class for inductor
class Inductor:
    # Constructor to intialize the object with input values
    def __init__(self, name, node1, node2, value):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.value = value

    # Function to modify the conducatnce matrix according to inductor in Modified nodal analysis
    def fillMatrix(self, G_matrix, nodes):
        k = nodes[self.node1]
        l = nodes[self.node2]

        if self.node1 == "GND":
            G_matrix[l][l] += 1.0 / self.value

        elif self.node2 == "GND":
            G_matrix[k][k] += 1.0 / self.value

        else:
            G_matrix[k][k] += 1.0 / self.value
            G_matrix[k][l] += -1.0 / self.value
            G_matrix[l][k] += -1.0 / self.value
            G_matrix[l][l] += 1.0 / self.value


# Class for Volatge Source
class VoltageSource:
    # Constructor to intialize the object with input values
    def __init__(self, name, node1, node2, value, count):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.value = value
        self.index = (
            count - 1
        )  # To store index of that Volatge source for future use in conductance matrix

    # Function to modify the conducatnce matrix and I vector according to Voltage source in Modified nodal analysis
    def fillMatrix(self, G_matrix, I_vector, nodes):
        k = nodes[self.node1]
        l = nodes[self.node2]
        vkl = node_count - 1 + self.index

        if self.node1 == "GND":
            G_matrix[vkl][l] -= 1
            G_matrix[l][vkl] -= 1

        elif self.node2 == "GND":
            G_matrix[vkl][k] -= 1
            G_matrix[k][vkl] -= 1

        else:
            G_matrix[k][vkl] += 1
            G_matrix[l][vkl] -= 1
            G_matrix[vkl][k] += 1
            G_matrix[vkl][l] -= 1

        I_vector[vkl][0] += self.value


# Independent Current source
class CurrentSource:
    # Constructor to intialize the object with input values
    def __init__(self, name, node1, node2, value):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.value = value

    # Function to modify the I vector according to Current Source in Modified nodal analysis
    def fillMatrix(self, I_vector, nodes):
        l = nodes[self.node1]
        k = nodes[self.node2]
        if self.node1 == "GND":
            I_vector[k] += self.value
        elif self.node2 == "GND":
            I_vector[l] += -self.value
        else:
            I_vector[l] += -self.value
            I_vector[k] += +self.value


# Voltage controlled voltage source
class VCVS:
    # node1, node2 = output
    # nc1, nc2 = control
    def __init__(self, name, node1, node2, nc1, nc2, gain, index):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.nc1 = nc1
        self.nc2 = nc2
        self.gain = gain
        self.index = index - 1

    # Update conducatnce matrix
    def fillMatrix(self, G_matrix, nodes):
        k = nodes[self.node1]
        l = nodes[self.node2]
        m = nodes[self.nc1]
        n = nodes[self.nc2]
        vkl = node_count - 1 + self.index

        if k != "GND":
            G_matrix[k][vkl] += 1
            G_matrix[vkl][k] += 1

        if l != "GND":
            G_matrix[l][vkl] += -1
            G_matrix[vkl][l] += -1

        if m != "GND":
            G_matrix[vkl][m] += -self.gain

        if n != "GND":
            G_matrix[vkl][n] += self.gain


# Voltage controlled current source
class VCCS:
    def __init__(self, name, node1, node2, nc1, nc2, gain):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.nc1 = nc1
        self.nc2 = nc2
        self.gain = gain

    # Update conductance matrix
    def fillMatrix(self, G_matrix, nodes):
        k = nodes[self.node1]
        l = nodes[self.node2]
        m = nodes[self.nc1]
        n = nodes[self.nc2]

        if k != "GND" and m != "GND":
            G_matrix[k][m] += self.gain

        if k != "GND" and n != "GND":
            G_matrix[k][n] += -self.gain

        if l != "GND" and m != "GND":
            G_matrix[l][m] += -self.gain

        if l != "GND" and n != "GND":
            G_matrix[l][n] += self.gain


class CCCS:
    def __init__(self, name, node1, node2, vcontrol, gain, index, values):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.vcontrol = vcontrol
        self.gain = gain
        self.index = index - 1
        self.nc1 = values[1]
        self.nc2 = values[2]


    def fillMatrix(self, G_matrix, nodes, values):
        nc1 = values[1]
        nc2 = values[2]

        k = nodes[self.node1]
        l = nodes[self.node2]
        m = nodes[nc1]
        n = nodes[nc2]
        imn = node_count - 1 + VS_count + self.index

        if k != "GND":
            G_matrix[k][imn] += self.gain

        if l != "GND":
            G_matrix[l][imn] -= self.gain

        if m != "GND":
            G_matrix[m][imn] += 1
            G_matrix[imn][m] += 1

        if n != "GND":
            G_matrix[n][imn] -= 1
            G_matrix[imn][n] -= 1


class CCVS:
    def __init__(self, name, node1, node2, vcontrol, gain, index, values):
        self.name = name
        self.node1 = node1
        self.node2 = node2
        self.vcontrol = vcontrol
        self.gain = gain
        self.index = index - 1
        self.nc1 = values[1]
        self.nc2 = values[2]

    def fillMatrix(self, G_matrix, nodes, values):
        nc1 = values[1]
        nc2 = values[2]
        k = nodes[self.node1]
        l = nodes[self.node2]
        m = nodes[nc1]
        n = nodes[nc2]
        imn = node_count - 1 + VS_count + self.index
        ikl = node_count - 1 + VS_count + self.index + 1

        if k != "GND":
            G_matrix[ikl][k] += 1
            G_matrix[k][ikl] += 1

        if m != "GND":
            G_matrix[imn][m] += 1
            G_matrix[m][imn] += 1

        if n != "GND":
            G_matrix[imn][n] -= 1
            G_matrix[n][imn] -= 1

        if l != "GND":
            G_matrix[ikl][l] -= 1
            G_matrix[l][ikl] -= 1

        G_matrix[ikl][imn] -= self.gain


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

    return words


# Store each component in their Class with appropriate value
def getComponents(LinesTokens, VI_count):
    components = []
    IC_count = 1
    VS_count = 1

    for line in LinesTokens:
        object_name = line[0]
        object = []
        name, n1, n2, *additional = line

        if object_name[0] == "R":  # Resistor
            value = additional[0]
            value = complex(float(value), 0)
            object = Resistor(name, n1, n2, value)

        elif object_name[0] == "C":  # Capacitor
            if ac_simulation == True:
                value = additional[0]
                value = complex(0, -1.0 / (omega * float(value)))
                object = Capacitor(name, n1, n2, value)

        elif object_name[0] == "L":  # Inductor
            value = additional[0]
            value = complex(0, omega * float(value))
            object = Inductor(name, n1, n2, value)

        elif object_name[0] == "V":  # Voltage Source
            sim = additional[0]
            value = 0
            if sim.upper() == "DC":  # Only amplitude
                amplitude = float(additional[1])
                value = amplitude
            elif sim.upper() == "AC":  # amplitude + phase
                amplitude = float(additional[1])
                phase = radians(float(additional[2]))
                value = (amplitude / 2) * complex(cos(phase), sin(phase))
            else:
                value = float(additional[0])
            object = VoltageSource(name, n1, n2, value, VS_count)
            VS_count += 1

        elif object_name[0] == "I":  # Current Source
            sim = additional[0]
            value = 0
            if sim.upper() == "DC":  # amplitude
                amplitude = float(additional[1])
                value = amplitude
            elif sim.upper() == "AC":  # amplitude + phase
                amplitude = float(additional[1])
                phase = radians(float(additional[2]))
                value = (amplitude / 2) * complex(cos(phase), sin(phase))
            object = CurrentSource(name, n1, n2, value)

        elif object_name[0] == "E":  # Voltage controlled Voltage source
            name, node1, node2, nc1, nc2, gain = line
            object = VCVS(name, node1, node2, nc1, nc2, float(gain), VS_count)
            VS_count += 1

        elif object_name[0] == "G":  # Voltage controlled current source
            name, node1, node2, nc1, nc2, gain = line
            object = VCCS(name, node1, node2, nc1, nc2, float(gain))

        elif object_name[0] == "F":  # Current controlled current source
            name, node1, node2, vcontrol, gain = line
            for element in components:
                if element.name == vcontrol:
                    values = [element.index, element.node1, element.node2]
            object = CCCS(name, node1, node2, vcontrol, float(gain), IC_count, values)
            IC_count += 1

        elif object_name[0] == "H":  # Current controlled voltage source
            name, node1, node2, vcontrol, gain = line
            for element in components:
                if element.name == vcontrol:
                    values = [element.index, element.node1, element.node2]
            object = CCCS(name, node1, node2, vcontrol, float(gain), IC_count, values)
            IC_count += 2

        components.append(object)

    VI_count[0] = VS_count - 1
    VI_count[1] = IC_count - 1

    return components


def addNode(nodes, count, node):
    if node not in nodes:
        if node == "GND":
            nodes[node] = "GND"
        else:
            nodes[node] = count
            count += 1

    return count


# Assign indices to node names
def assignNodeNames(components):
    global node_count
    nodes = {}
    count = 0
    for object in components:
        name = object.name

        # Add nodes to dictionary
        node1, node2 = object.node1, object.node2
        count = addNode(nodes, count, node1)
        count = addNode(nodes, count, node2)

        # Extra node inputs in line if it is a controlled source
        if name[0] == "E" or name[0] == "G" or name[0] == "F" or name[0] == "H":
            nc1, nc2 = object.nc1, object.nc2
            count = addNode(nodes, count, nc1)
            count = addNode(nodes, count, nc2)

    node_count = len(nodes)

    return nodes


# Find number of Voltage Sources and VCVS in Input
def getVSCount(components):
    cnt = 0
    for object in components:
        name = type(object).__name__
        if name == V or name == E:
            cnt = cnt + 1

    return cnt


# To Update Conductance matrix
def updateMatrix(G_matrix, I_vector, components, nodes):
    for object in components:
        if type(object).__name__ == R:
            object.fillMatrix(G_matrix, nodes)

        elif type(object).__name__ == V:
            object.fillMatrix(G_matrix, I_vector, nodes)

        elif type(object).__name__ == I:
            object.fillMatrix(I_vector, nodes)

        elif type(object).__name__ == C:
            object.fillMatrix(G_matrix, nodes)

        elif type(object).__name__ == L:
            object.fillMatrix(G_matrix, nodes)

        elif type(object).__name__ == E:
            object.fillMatrix(G_matrix, nodes)

        elif type(object).__name__ == G:
            object.fillMatrix(G_matrix, nodes)

        elif type(object).__name__ == F:
            for element in components:
                if element.name == object.vcontrol:
                    values = [element.index, element.node1, element.node2]
            object.fillMatrix(G_matrix, nodes, values)

        elif type(object).__name__ == H:
            for element in components:
                if element.name == object.vcontrol:
                    values = [element.index, element.node1, element.node2]
            object.fillMatrix(G_matrix, nodes, values)


# Print Volatge and Current outputs
def printResult(V_vector, nodes, components):
    print("Reference: GND = 0 V", end="\n\n")
    for node in nodes:
        if node != "GND":
            print(f"V at node {node}: {V_vector[nodes[node]]}", end="\n\n")

    for object in components:
        if type(object).__name__ == V:
            print(
                f"I through {object.name}: {V_vector[node_count - 1 + object.index]}",
                end="\n\n",
            )

        if type(object).__name__ == E:
            print(
                f"Current through {object.name}: {V_vector[node_count - 1 + object.index]}",
                end="\n\n",
            )

        if type(object).__name__ == F:
            print(
                f"Current through {object.name}: {V_vector[node_count - 1 + VS_count + object.index]}",
                end="\n\n",
            )


if __name__ == "__main__":

    if len(argv) != 2:  # If more than 2 arguments are passed
        print("Invalid operation !")
        print(f"Usage: {argv[0]} <inputfile>'")
        exit()

    try:
        with open(argv[1]) as f:
            lines = f.readlines()
            start = -1
            end = -2
            ac = -1
            start_cnt = 0
            end_cnt = 0
            for line in lines:  # extracting circuit definition: start and end lines
                if CIRCUIT == line[0 : len(CIRCUIT)]:
                    start = lines.index(line)
                    start_cnt = start_cnt + 1
                elif END == line[: len(END)]:
                    end = lines.index(line)
                    end_cnt = end_cnt + 1
                elif AC == line[: len(AC)]:
                    ac = lines.index(line)

            if (
                start_cnt != 1 or end_cnt != 1
            ):  # Check if there are multiple .circuit/.end declarations in input file
                print(
                    "Invalid circuit definition!\nMultiple/No .circuit/.end declarations detected"
                )
                exit()

            if start >= end:  # validating circuit block
                print("Invalid circuit definition!")
                exit()

            if ac != -1:
                ac_simulation = True
                _, _, frequency = lines[ac].split("#")[0].split()
                omega = 2.0 * PI * float(frequency)

            LinesRaw = lines[start + 1 : end]
            LinesActual = extractData(LinesRaw)

            # Split values in each line
            LinesTokens = [extractTokens(line) for line in LinesActual]

            # Get list of objects with appropriate classes
            components = getComponents(LinesTokens, VI_count)
            VS_count = VI_count[0]
            IC_count = VI_count[1]

            # Get dictionary of assigned node names
            
            nodes = assignNodeNames(components)
            node_count = len(nodes)

            # Get Voltage sources count

            rows = node_count +  VS_count + IC_count - 1
            columns = rows

            # Resize Matrices according to input
            G_matrix.resize((rows, columns))
            I_vector.resize((columns, 1))

            updateMatrix(G_matrix, I_vector, components, nodes)

            print(G_matrix)
            V_vector = np.linalg.solve(G_matrix, I_vector)

            
            print(I_vector)
            printResult(V_vector, nodes, components)

    except IOError:
        print("Invalid file!")
        exit()
