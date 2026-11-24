
import re
from openpyxl import Workbook

# Input and output file paths
from pathlib import Path

# folder of the current Python script
cartella = Path(__file__).resolve().parent

# Input and output file paths, in the same folkder of the script
input_file = cartella / "report.txt"
output_file = cartella / "results.xlsx"

# the 25 metrics to extract from the report
metriche = [
    "Signal-to-Interference-plus-Noise Ratio",
    "Error Vector Magnitude percentage",
    "Error Vector Magnitude decibels",
    "E_b/N_0 linear",
    "E_b/N_0 decibels",
    "Hard errors inferred from LLR sign",
    "Total final errored bits from original packets",
    "Bit Error Rate before decoding",
    "Bit Error Rate after decoding",
    "Total bursts erased",
    "Total packets completely reconstructed",
    "Burst Loss Rate",
    "Packet Error Rate",
    "Total number of detections made by CA-CFAR",
    #"Total number of detections made by FIXED DETECTION THRESHOLD",
    "Total TRUE detections",
    "True Transmitter detections",
    "True Friends (Interferers) detections",
    "Total FALSE detections",
    "False Detection Probability (FDP)",
    "True Positive windows",
    "False Positive windows",
    "False Negative windows",
    "True Negative windows",
    "Precision (TP/(TP+FP))",
    "Recall (TP/(TP+FN))"
]

# Read the report
with open(input_file, "r", encoding="utf-8") as f:
    righe = f.readlines()

valori = []

# Extracts the first number of each metric from the report
for metrica in metriche:
    trovato = False

    for riga in righe:
        if riga.strip().startswith(metrica):
            parte_numerica = riga.split("=", 1)[1]

            match = re.search(
                r"[-+]?\d+(?:\.\d+)?",
                parte_numerica
            )

            if match:
                numero = float(match.group())
                valori.append(numero)
                trovato = True
                break

    if not trovato:
        raise ValueError(f"Metric not found: {metrica}")

# Create the Excel file
wb = Workbook()
ws = wb.active
ws.title = "Results"

# Write the 25 numbers horizontally in the first row
for colonna, valore in enumerate(valori, start=1):
    ws.cell(row=1, column=colonna, value=valore)

# Save the file Excel
wb.save(output_file)

print(f"Created file: {output_file}")
print(f"Extracted values: {len(valori)}")