#!/usr/bin/env python3
"""Parse the processed CSVs and regenerate selected plots.
Run this from the repository root after installing pandas and matplotlib.
"""
# The full plotting script used to generate the distributed figures is in the ChatGPT artifact history.
# This lightweight placeholder documents the input files and can be expanded as needed.
from pathlib import Path
import pandas as pd

base = Path('data/processed')
print(pd.read_csv(base / 'summary_percentage.csv').head())
