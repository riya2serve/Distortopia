#!/usr/bin/env python3
"""Main CLI interface for Distortopia."""

import argparse
from distortopia import simf1poly, segdistorters

def main():
    parser = argparse.ArgumentParser(description="Distortopia CLI")
    parser.add_argument("--simulate-f1", action="store_true", help="Simulate F1 hybrid")
    parser.add_argument("--detect-distortion", action="store_true", help="Detect segregation distortion")
    # More arguments here...

    args = parser.parse_args()

    if args.simulate_f1:
        __simf1poly__.run_simulation(args) #will either simulate hybrids
    elif args.detect_distortion: #otherwise will detect distortions
        __segdistorters__.detect(args)

if __name__ == "__main__":
    main()
