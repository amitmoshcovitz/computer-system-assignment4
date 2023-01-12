import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Optional

DEFAULT_RUN_COUNT = 3

TEST_SERVER = "planet"

RUN_COUNT = int(sys.argv[1]) if len(sys.argv) == 2 else DEFAULT_RUN_COUNT
HOST = subprocess.check_output("hostname").decode().strip()

SLEEP_PERIOD = 0.4 if os.system("timeout 0.25 ls > /dev/null 2>/dev/null") == 0 else 1

TEMP_FILE_PATH = Path(".compare_tmp")

ORIGINAL_FILES = [
    Path("200_pic/Heritage_200.bmp"),
    Path("500_pic/gibson_500.bmp"),
    Path("700_pic/Matera_700.bmp")
]

PARAMS = ["1", "2", "3"]
PARAMS_TO_FORM_FIELDS = {
    "1": "entry.295728809",
    "2": "entry.266687381",
    "3": "entry.296772969"
}

RESULT_TYPES = ["Blur", "Filtered_Blur", "Filtered_Sharpen", "Row_Blur", "Row_Sharpen", "Sharpen"]

RESULTS = dict({
    str(filename): dict({param: [] for param in PARAMS}) for filename in ORIGINAL_FILES
})

def print_with_emoji(msg):
    sys.stdout.buffer.write(f"{msg}\n".encode())

def run_program(full_path, param) -> Optional[float]:
    print(f"Running: {full_path} {param} ... ", end='', flush=True)
    os.system(f'timeout {SLEEP_PERIOD} stdbuf -oL ./showBMP "{full_path}" {param} > {TEMP_FILE_PATH}')
    output = TEMP_FILE_PATH.read_text()
    if TEMP_FILE_PATH.exists():
        TEMP_FILE_PATH.unlink()
    try:
        return float(output.splitlines()[-1].split()[2])  # Message is of format: Total runtime: 5.128000 ms
    except:
        # In case of an error, just return None
        return None


def clear_results():
    for result_type in RESULT_TYPES:
        os.system(f"rm -f {result_type}.bmp")


def main():
    ERROR = False
    print("Rebuilding...")
    os.system("make clean")
    os.system("make")
    print()

    print(f"Running test {RUN_COUNT} times...")

    for iteration in range(RUN_COUNT):
        print("*" * 100)
        print(f"    ITERATION {iteration + 1}    ".center(100, "*"))
        print("*" * 100, end="\n\n")
        for full_path in ORIGINAL_FILES:
            path_str = str(full_path)
            clear_results()

            for param in PARAMS:
                time_took = run_program(path_str, param)
                if time_took is None:
                    print_with_emoji("❌ Failed (Error)")
                    ERROR = True
                    continue
                else:
                    print(f"Took {time_took:.3f} ms")
                    RESULTS[path_str][param].append(time_took)
                    print()

            for result_type in RESULT_TYPES:
                print(f"Checking {result_type} ... ".ljust(50), end='', flush=True)
                correct = os.system(f"cmp {result_type}.bmp {str(full_path.parent)}/{result_type}_correct.bmp > /dev/null 2>&1") == 0
                if correct:
                    print_with_emoji("✅ Passed")
                else:
                    print_with_emoji("❌ Failed (Wrong output)")
                    ERROR = True
                print()

    print("  Runtimes  ".center(100, "="))

    for full_path, timings in RESULTS.items():
        print(f"    {full_path:20s}    ".center(100, "="))
        for param, timing in timings.items():
            average = sum(timing) / RUN_COUNT
            print(f"Param: {param} Average runtime: {average:.3f}")

    print()

    if HOST != TEST_SERVER:
        print_with_emoji(f"❗ Not uploading results because not running on {TEST_SERVER}")
    elif ERROR:
        print_with_emoji("❗ Not uploading results because not all runs were successful")
    elif RUN_COUNT != 20:
        print_with_emoji("❗ Not uploading results because the average wasn't on 20 runs")
    else:
        print_with_emoji("📤 Uploading results to Google Forms...")

        entries_data = []
        for param in PARAMS:
            results = []
            for filename in ORIGINAL_FILES:
                results += RESULTS[str(filename)][param]
            average = sum(results) / float(len(results))

            entries_data.append(f"{PARAMS_TO_FORM_FIELDS[param]}={average:.3f}")

        os.system(f"""
    curl 'https://docs.google.com/forms/d/e/1FAIpQLScIN4njh40n6osL4Ag4qJf5RsS88_q6NqEHYwWJeI-Wn1seJw/formResponse' \
      -H 'authority: docs.google.com' \
      -H 'accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9' \
      -H 'accept-language: en-GB,en;q=0.9' \
      -H 'cache-control: max-age=0' \
      -H 'content-type: application/x-www-form-urlencoded' \
      -H 'origin: https://docs.google.com' \
      -H 'referer: https://docs.google.com/forms/d/e/1FAIpQLScIN4njh40n6osL4Ag4qJf5RsS88_q6NqEHYwWJeI-Wn1seJw/viewform?fbzx=4387170321690009457' \
      -H 'sec-ch-ua: ".Not/A)Brand";v="99", "Google Chrome";v="103", "Chromium";v="103"' \
      -H 'sec-ch-ua-arch: "x86"' \
      -H 'sec-ch-ua-bitness: "64"' \
      -H 'sec-ch-ua-full-version-list: ".Not/A)Brand";v="99.0.0.0", "Google Chrome";v="103.0.5060.53", "Chromium";v="103.0.5060.53"' \
      -H 'sec-ch-ua-mobile: ?0' \
      -H 'sec-ch-ua-model: ""' \
      -H 'sec-ch-ua-platform: "Linux"' \
      -H 'sec-ch-ua-platform-version: "5.15.0"' \
      -H 'sec-ch-ua-wow64: ?0' \
      -H 'sec-fetch-dest: document' \
      -H 'sec-fetch-mode: navigate' \
      -H 'sec-fetch-site: same-origin' \
      -H 'sec-fetch-user: ?1' \
      -H 'upgrade-insecure-requests: 1' \
      -H 'user-agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/103.0.0.0 Safari/537.36' \
      --data "{'&'.join(entries_data)}&fvv=1&partialResponse=%5Bnull%2Cnull%2C%224387170321690009457%22%5D&pageHistory=0&fbzx=4387170321690009457" \
          --compressed > /dev/null
        """)


    print("View others' results at: https://docs.google.com/spreadsheets/d/1IRb06fPS6NSGQOE4hzGyIVe3N11uUlqWLj6H94guaZY")
    print("Test script was written by Ron Even")
    print("Thank you for using it!")

main()

