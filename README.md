## Para primeira execução:

``
pip install -r requirements.txt
``

``
pip install -e .
``

## Run via código:
```python
from run.Simulation import Simulation

Simulation.startSimulation(mode='vad', inputPath='./input.txt')

````

## Run via terminal:
``
python init.py -m vad .\input.txt
``
