from qiskit.circuit import Gate
from qiskit.circuit import QuantumCircuit
from qiskit.circuit import QuantumRegister
from qiskit.extensions.standard import HGate, Cu1Gate, U1Gate, CnotGate, XGate
from qiskit.extensions.standard import FredkinGate

import numpy as np

##################################### QFT #####################################

class QFTGate(Gate):
    """QFT gate."""

    def __init__(self, n, inverse=False):
        """Create new QFT gate."""
        if inverse:
            super().__init__("iQFT", n, [])
        else:
            super().__init__("QFT", n, [])
        self.n = n
        self._inverse = inverse

    def _define(self):
        definition = []
        q = QuantumRegister(self.n, "q")
        rule = []

        if not self._inverse:
            for j in range(self.n - 1, -1, -1):
                rule.append((HGate(), [q[j]], []))
                for i in range(j):
                    rule.append((Cu1Gate(np.pi/float(2**(i + 1))), [q[j - (i + 1)], q[j]], []))
        else:
            for j in range(self.n):
                for i in range(j - 1, -1, -1):
                    rule.append((Cu1Gate(-np.pi/float(2**(i + 1))), [q[j - (i + 1)], q[j]], []))
                rule.append((HGate(), [q[j]], []))

        for inst in rule:
            definition.append(inst)
        self.definition = definition

    def inverse(self):
        """Invert this gate."""
        return QFTGate(self.n, inverse=not self._inverse)

def qft(self, tgt, n):
    """Apply QFT to tgt."""
    return self.append(QFTGate(n), [tgt[i] for i in range(n)], [])

def iqft(self, tgt, n):
    """Apply iQFT to tgt."""
    return self.append(QFTGate(n, inverse=True), [tgt[i] for i in range(n)], [])


QuantumCircuit.qft = qft
QuantumCircuit.iqft = iqft


################################### PhiADDa ###################################
ADD_STRING = """Phi
ADD
(%d)"""
SUB_STRING = """Phi
SUB
(%d)"""

class PhiADDaGate(Gate):
    """PhiADDa gate."""

    def __init__(self, a, n, c1=None, c2=None, inverse=False):
        """Create new PhiADDa gate."""
        self.total_n = n
        if c1:
            self.total_n += 1
        if c2:
            self.total_n += 1

        if inverse:
            super().__init__(SUB_STRING % a, self.total_n, [])
        else:
            super().__init__(ADD_STRING % a, self.total_n, [])
        self.n = n
        self.a = a
        self.c1 = c1
        self.c2 = c2
        self._inverse = inverse

    def _define(self):
        definition = []
        # No controls:
        #   q = b
        # One Control:
        #   q = b, c1
        # Two Controls:
        #   q = c1, c2, b
        q = QuantumRegister(self.total_n, "q")
        rule = []

        angles = [0]*self.n
        for i in range(self.n):
            if self.a & (1 << i):
                for j in range(i, self.n):
                    angles[j] += np.pi/2**(j-i)

        # Inverse of U1 gates is just a negative theta
        if self._inverse:
            for i in range(len(angles)):
                angles[i] *= -1

        # No controlled bits
        if self.c1 is None and self.c2 is None:
            for i in range(len(angles)):
                rule.append((U1Gate(angles[i]), [q[i]], []))

        # One controlled bit
        if self.c1 and self.c2 is None:
            for i in range(len(angles)):
                rule.append((Cu1Gate(angles[i]), [q[self.total_n - 1], q[i]], []))

        # Two controlled bits
        # Uses sqrt of U1 gates which is just half of theta
        if self.c1 and self.c2:
            for i in range(len(angles)):
                rule.append((Cu1Gate(angles[i]/2.), [q[1], q[i + 2]], []))
                # circ.cu1(angles[i]/2., c2, b[i])
            rule.append((CnotGate(), [q[0], q[1]], []))
            # circ.cx(c1, c2)

            for i in range(len(angles)):
                rule.append((Cu1Gate(-angles[i]/2.), [q[1], q[i + 2]], []))
                # circ.cu1(-angles[i]/2., c2, b[i])
            rule.append((CnotGate(), [q[0], q[1]], []))
            # circ.cx(c1, c2)
            for i in range(len(angles)):
                rule.append((Cu1Gate(angles[i]/2.), [q[0], q[i + 2]], []))
                # circ.cu1(angles[i]/2., c1, b[i])

        for inst in rule:
            definition.append(inst)
        self.definition = definition

    def inverse(self):
        """Invert this gate."""
        return PhiADDaGate(self.a, self.n, c1=self.c1, c2=self.c2, inverse=not self._inverse)

def PhiADDa(self, a, b, n):
    """Apply a plus b."""
    return self.append(PhiADDaGate(a, n), [b[i] for i in range(n)], [])

def cPhiADDa(self, a, b, c1, n):
    """Apply a plus b if c1==1."""
    return self.append(PhiADDaGate(a, n, c1=c1), [b[i] for i in range(n)] + [c1[0]], [])

def ccPhiADDa(self, a, b, c1, c2, n):
    """Apply a plus b if c1==c2==1."""
    return self.append(PhiADDaGate(a, n, c1=c1, c2=c2), [c1, c2] + [b[i] for i in range(n)], [])

def PhiSUBa(self, a, b, n):
    """Apply b minus a."""
    return self.append(PhiADDaGate(a, n, inverse=True), [b[i] for i in range(n)], [])

def cPhiSUBa(self, a, b, c1, n):
    """Apply b minus a if c1==1."""
    return self.append(PhiADDaGate(a, n, c1=c1, inverse=True), [b[i] for i in range(n)] + [c1[0]], [])

def ccPhiSUBa(self, a, b, c1, c2, n):
    """Apply b minus a if c1==c2==1."""
    return self.append(PhiADDaGate(a, n, c1=c1, c2=c2, inverse=True), [c1, c2] + [b[i] for i in range(n)], [])

QuantumCircuit.PhiADDa = PhiADDa
QuantumCircuit.cPhiADDa = cPhiADDa
QuantumCircuit.ccPhiADDa = ccPhiADDa
QuantumCircuit.PhiSUBa = PhiSUBa
QuantumCircuit.cPhiSUBa = cPhiSUBa
QuantumCircuit.ccPhiSUBa = ccPhiSUBa

################################# PhiADDaModN #################################
ADD_MOD_STRING = """Phi
ADD
(%d)
Mod
%d"""
SUB_MOD_STRING = """Phi
SUB
(%d)
Mod
%d"""

class PhiADDaModNGate(Gate):
    """PhiADDaModN gate."""

    def __init__(self, a, N, n, inverse=False):
        """Create new PhiADDaModN gate."""
        self.total_n = n + 3

        if inverse:
            super().__init__(SUB_MOD_STRING % (a, N), self.total_n, [])
        else:
            super().__init__(ADD_MOD_STRING % (a, N), self.total_n, [])
        self.n = n
        self.N = N
        self.a = a
        self._inverse = inverse

    def _define(self):
        definition = []
        #   q = c1, c2, b, anc
        q = QuantumRegister(self.total_n, "q")

        # No easy way to inverse other than reverse the order of gates and
        # invert them.
        if not self._inverse:
            rule = [
                (PhiADDaGate(self.a, self.n, c1=q[0], c2=q[1]), [q[0], q[1]] + [q[i] for i in range(2, self.n + 2)], []),
                (PhiADDaGate(self.N, self.n, inverse=True), [q[i] for i in range(2, self.n + 2)], []),

                (QFTGate(self.n, inverse=True), [q[i] for i in range(2, self.n + 2)], []),
                (CnotGate(), [q[self.total_n - 2], q[self.total_n - 1]], []),
                (QFTGate(self.n), [q[i] for i in range(2, self.n + 2)], []),

                (PhiADDaGate(self.N, self.n, c1=q[self.total_n - 1]), [q[i] for i in range(2, self.n + 2)] + [q[self.total_n - 1]], []),

                (PhiADDaGate(self.a, self.n, c1=q[0], c2=q[1], inverse=True), [q[0], q[1]] + [q[i] for i in range(2, self.n + 2)], []),

                (QFTGate(self.n, inverse=True), [q[i] for i in range(2, self.n + 2)], []),

                (XGate(), [q[self.total_n - 2]], []),
                (CnotGate(), [q[self.total_n - 2], q[self.total_n - 1]], []),
                (XGate(), [q[self.total_n - 2]], []),

                (QFTGate(self.n), [q[i] for i in range(2, self.n + 2)], []),

                (PhiADDaGate(self.a, self.n, c1=q[0], c2=q[1]), [q[0], q[1]] + [q[i] for i in range(2, self.n + 2)], []),
            ]
        else:
            rule = [
                (PhiADDaGate(self.a, self.n, c1=q[0], c2=q[1], inverse=True), [q[0], q[1]] + [q[i] for i in range(2, self.n + 2)], []),

                (QFTGate(self.n, inverse=True), [q[i] for i in range(2, self.n + 2)], []),

                (XGate(), [q[self.total_n - 2]], []),
                (CnotGate(), [q[self.total_n - 2], q[self.total_n - 1]], []),
                (XGate(), [q[self.total_n - 2]], []),

                (QFTGate(self.n), [q[i] for i in range(2, self.n + 2)], []),

                (PhiADDaGate(self.a, self.n, c1=q[0], c2=q[1]), [q[0], q[1]] + [q[i] for i in range(2, self.n + 2)], []),

                (PhiADDaGate(self.N, self.n, c1=q[self.total_n - 1], inverse=True), [q[i] for i in range(2, self.n + 2)] + [q[self.total_n - 1]], []),

                (QFTGate(self.n, inverse=True), [q[i] for i in range(2, self.n + 2)], []),
                (CnotGate(), [q[self.total_n - 2], q[self.total_n - 1]], []),
                (QFTGate(self.n), [q[i] for i in range(2, self.n + 2)], []),

                (PhiADDaGate(self.N, self.n), [q[i] for i in range(2, self.n + 2)], []),
                (PhiADDaGate(self.a, self.n, c1=q[0], c2=q[1], inverse=True), [q[0], q[1]] + [q[i] for i in range(2, self.n + 2)], []),
            ]

        for inst in rule:
            definition.append(inst)
        self.definition = definition

    def inverse(self):
        """Invert this gate."""
        return PhiADDaModNGate(self.a, self.N, self.n, inverse=not self._inverse)

def PhiADDaModN(self, a, b, c1, c2, anc, N, n):
    """Apply (a + b) % N if c1==c2==1."""
    return self.append(PhiADDaModNGate(a, N, n), [c1, c2] + [b[i] for i in range(n)] + [anc], [])

def PhiSUBaModN(self, a, b, c1, c2, anc, N, n):
    """Apply (a - b) % N if c1==c2==1."""
    return self.append(PhiADDaModNGate(a, N, n, inverse=True), [c1, c2] + [b[i] for i in range(n)] + [anc], [])

QuantumCircuit.PhiADDaModN = PhiADDaModN
QuantumCircuit.PhiSUBaModN = PhiSUBaModN

################################# CMULTaModN ##################################
CMULT_STRING = """C
MULT
(%d)
Mod
%d"""
CDIV_STRING = """C
DIV
(%d)
Mod
%d"""

class CMULTaModNGate(Gate):
    """CMULTaModN gate."""

    def __init__(self, a, N, n, inverse=False):
        """Create new CMULTaModN gate."""
        self.total_n = 2*n + 2

        if inverse:
            super().__init__(CDIV_STRING % (a, N), self.total_n, [])
        else:
            super().__init__(CMULT_STRING % (a, N), self.total_n, [])
        self.n = n
        self.N = N
        self.a = a
        self._inverse = inverse

    def _define(self):
        definition = []
        #   q = c, x, b, anc
        q = QuantumRegister(self.total_n, "q")

        rule = [
            (QFTGate(self.n), [q[i] for i in range(self.n + 1, self.total_n - 1)], []),
        ]

        if self._inverse:
            bounds = range(self.n - 1, -1, -1)
        else:
            bounds = range(self.n)

        for i in bounds:
            rule.append((PhiADDaModNGate(((2**i)*self.a) % self.N, self.N, self.n, inverse=self._inverse), [q[0], q[i + 1]] + [q[j] for j in range(self.n + 1, self.total_n - 1)] + [q[self.total_n - 1]], []))


        rule.append((QFTGate(self.n, inverse=True), [q[i] for i in range(self.n + 1, self.total_n - 1)], []))

        for inst in rule:
            definition.append(inst)
        self.definition = definition

    def inverse(self):
        """Invert this gate."""
        return CMULTaModNGate(self.a, self.N, self.n, inverse=not self._inverse)

def CMULTaModN(self, a, b, c, x, anc, N, n):
    """Apply (b + a*x) % N if c==1."""
    return self.append(CMULTaModNGate(a, N, n), [c] + [x[i] for i in range(n)] + [b[i] for i in range(n)] + [anc], [])

def CDIVaModN(self, a, b, c, x, anc, N, n):
    """Apply (b - a*x) % N if c==1."""
    return self.append(CMULTaModNGate(a, N, n, inverse=True), [c] + [x[i] for i in range(n)] + [b[i] for i in range(n)] + [anc], [])


QuantumCircuit.CMULTaModN = CMULTaModN
QuantumCircuit.CDIVaModN = CDIVaModN

#################################### Cua ######################################
def modInverse(a, m):
    """
    Finds the modular multiplicative inverse of a mod m.
    a and m must be coprime.

    modInverse of a is x such that a*x == 1 (mod m)
    Ex. a = 27, m=392
    x = 363 since 27*363 == 1 (mod 392)
    """
    m0 = m
    y = 0
    x = 1

    if (m == 1) :
        return 0

    while (a > 1) :
        # q is quotient
        q = a // m

        t = m

        # m is remainder now, process
        # same as Euclid's algo
        m = a % m
        a = t
        t = y

        # Update x and y
        y = x - q * y
        x = t

    # Make x positive
    if (x < 0) :
        x = x + m0

    return x

class CuaGate(Gate):
    """Cua gate."""

    def __init__(self, a, N, n):
        """Create new Cua gate."""
        self.total_n = 2*n + 2
        super().__init__("U %d" % (a,), self.total_n, [])
        self.n = n
        self.N = N
        self.a = a

    def _define(self):
        definition = []
        #   q = c, x, z, anc
        q = QuantumRegister(self.total_n, "q")

        rule = [
            (CMULTaModNGate(self.a, self.N, self.n), [q[i] for i in range(self.total_n)], []),
        ]

        for i in range(self.n):
            rule.append((FredkinGate(), [q[0], q[i + 1], q[self.n + 1 + i]], [])),

        rule.append((CMULTaModNGate(modInverse(self.a, self.N), self.N, self.n, inverse=True), [q[i] for i in range(self.total_n)], []))


        for inst in rule:
            definition.append(inst)
        self.definition = definition

def cua(self, a, c, x, z, anc, N, n):
    """Apply (a*x) % N if c==1."""
    return self.append(CuaGate(a, N, n), [c] + [x[i] for i in range(n)] + [z[i] for i in range(n)] + [anc], [])

QuantumCircuit.cua = cua
