#! /bin/env python

# Converts PROSITE pattern to Regular Expression

def convertPattern(pattern):
	# Replace " " with nothing
    pattern = pattern.replace(" ","")
    # Replace "-" with nothing
    pattern = pattern.replace("-","")
    # Replace last "." with nothing
    pattern = pattern.replace(".","")
    # Replace "<" with "^"
    pattern = pattern.replace("<","^")
    # Replace ">" with "$"
    pattern = pattern.replace(">","$")
    # Replace "{}" with "[^]"
    pattern = pattern.replace("{","[^")
    pattern = pattern.replace("}","]")
    # Replace "()" with "{}"
    pattern = pattern.replace("(","{")
    pattern = pattern.replace(")","}")
    # Replace "x" with "."
    pattern = pattern.replace("x",".")
    # All caps
    pattern = pattern.upper()
    return pattern

if __name__ == "__main__":
    import sys

    if len(sys.argv[1:]):
        for line in sys.argv[1:]:
            print convertPattern(line.strip())
    else:
        line = sys.stdin.readline()
        while line:
            print convertPattern(line.strip())
            line = sys.stdin.readline()
