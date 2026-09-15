#!/usr/bin/env python3
import os
import sys

if __package__ is None or __package__ == "":
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from conservation.postprocess import main


if __name__ == "__main__":
    sys.argv.insert(1, "--mode")
    sys.argv.insert(2, "pc")
    main()

