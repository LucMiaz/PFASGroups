#!/usr/bin/env python3

# Script to fix the generate_unified_report.py by removing duplicate sections

with open('generate_unified_report.py', 'r') as f:
    lines = f.readlines()

# Find the first main function
first_main_idx = None
second_main_idx = None

for i, line in enumerate(lines):
    if 'def main():' in line:
        if first_main_idx is None:
            first_main_idx = i
        else:
            second_main_idx = i
            break

print(f"First main function at line {first_main_idx + 1}")
print(f"Second main function at line {second_main_idx + 1}")

# Keep everything up to the first main function (exclusive) and from the second main function onwards
fixed_lines = lines[:first_main_idx] + lines[second_main_idx:]

# Write the fixed content
with open('generate_unified_report.py', 'w') as f:
    f.writelines(fixed_lines)

print(f"Fixed file. Removed {len(lines) - len(fixed_lines)} lines.")
print(f"Original: {len(lines)} lines, Fixed: {len(fixed_lines)} lines")