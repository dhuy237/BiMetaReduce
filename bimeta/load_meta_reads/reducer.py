#!/usr/bin/env python3

from operator import itemgetter 
import sys

outKey = None # output Key   (print to STDOUT)   
outValue = None # output Value (print to STDOUT)
key = None
valueNumber = None
valueRead = None

for line in sys.stdin:
    line = line.strip() # remove the whitespace at the beginning and at the end
    key, value = line.split() # split at space
    
    # Check if "value" variable is label (a number) or a read; by converting value (currently a string) to int
    try: 
        value = int(value)
        valueNumber = value
    except ValueError:  # is not a number;
        valueRead = value
        #continue
        
    # Hadoop shuffle/sort phase in Reducer sort by the Key
    if outKey == key:
        outValue = str(valueNumber) + " " + valueRead    
    else:
        if outKey:
            # Write result to Standard Ouput STDOUT
            print('%s %s' % (outKey, outValue))
        outKey = key
        
if outKey == key:
    print('%s %s' % (outKey, outValue))