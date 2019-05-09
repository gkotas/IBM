# IBM Qiskit

## Setup the Environment
Python 3.6.7

Create a virtualenv
```
virtualenv env
```

Source the environment (in Linux)
```
source env/bin/activate
```

Install Qiskit and the dependances
```
pip install -r requirements.txt
```

Edit Qconfig.py and add your own API Token.

## Shor's in 2n+3
An implementation of Shor's algorithm in just 2n+3 qubits. From this paper https://arxiv.org/abs/quant-ph/0205095v3 by Stephane Beauregard.

### Custom Gates
Uses custom gates for nicer looking circuit drawings. Builds up all of the necessary gates to create a Controlled U-a gate.

### No Custom Gates
Just uses functions to create all of the gates using the basic qiskit gates. Same functionality as the Custom Gate implementations but without the nice circuit drawings.

## Quantum Adder
A notebook to demonstrate the differences between a simulation, noisy simulation, and a real run. Uses the quantum addition circuit from the Custom Gates.

