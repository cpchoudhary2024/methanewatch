"""Subset the Flaring Monitor observation table to the Permian box."""

import sys

import pandas as pd

LON0, LAT0, LON1, LAT1 = -105.0, 30.8, -101.0, 33.5

src, dst = sys.argv[1:3]
df = pd.read_csv(src, low_memory=False)
keep = df.longitude.between(LON0, LON1) & df.latitude.between(LAT0, LAT1)
df[keep].to_csv(dst, index=False)
