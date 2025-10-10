import os
import sys
import requests

def download_file(url, local_filename):
    """Download a file from a URL to a local path."""
    try:
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(local_filename, 'wb') as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        return local_filename
    except Exception as e:
        print(f"Error downloading file: {e}")
        return None

def main():
    print(f"Current working directory: {os.getcwd()}")

    input_path = input("Enter input file path/filename or URL: ").strip()

    # Check if input is URL (very simple check)

    try:
        with open('s1223.dat', 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) < 2:
            print("Input file doesn't have enough lines to process (needs at least 2).")
            return

        # Skip first line, process the rest
        data_lines = lines[1:]

        processed = []
        for line in data_lines:
            parts = line.split()
            if len(parts) == 2:
                parts.append('0')
                processed.append(parts)
            else:
                print(f"Warning: skipping malformed line: {line}")


    except Exception as e:
        print(f"An error occurred during processing: {e}")

if __name__ == "__main__":
    main()
