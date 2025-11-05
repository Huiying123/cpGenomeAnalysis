import os
from Bio import SeqIO

# --- User Configuration ---
GENBANK_DIR = "/data01/wangkun/projects/taibaishaolan/03.analysis/genbank/" #<-- EDIT THIS
OUTPUT_FILE = "junction_summary.tsv"
# -------------------------

def is_inverted_repeat(feature):
    """Checks if a feature is an Inverted Repeat based on several common annotation styles."""
    if feature.type != "repeat_region":
        return False
    if "rpt_family" in feature.qualifiers and feature.qualifiers.get('rpt_family') == ['Inverted Repeat']:
        return True
    if "note" in feature.qualifiers and "inverted repeat" in feature.qualifiers.get('note')[0].lower():
        return True
    if "rpt_type" in feature.qualifiers and feature.qualifiers.get('rpt_type') == ['inverted']:
        return True
    return False

def get_plastome_structure(genbank_file):
    """Parses a GenBank file to calculate the lengths of LSC, SSC, and IR regions."""
    try:
        record = SeqIO.read(genbank_file, "genbank")
        total_length = len(record.seq)
        ir_regions = []

        # Find all annotated inverted repeat regions using the flexible checker
        for feature in record.features:
            if is_inverted_repeat(feature):
                ir_regions.append(feature.location)

        if len(ir_regions) < 2:
            return total_length, "N/A", "N/A", "N/A", "Error: Not enough IRs found"

        # Sort by start position to identify IRa and IRb
        ir_regions.sort(key=lambda loc: loc.start)
        ira = ir_regions[0]
        irb = ir_regions[1]

        # Calculate lengths
        ir_len = len(ira)
        lsc_len = irb.start - ira.end
        ssc_len = total_length - lsc_len - (2 * ir_len)
        
        return total_length, lsc_len, ssc_len, ir_len, "Success"

    except Exception as e:
        return "ERROR", "ERROR", "ERROR", "ERROR", f"Error: {e}"

# --- Main Script Execution ---
if __name__ == "__main__":
    files_to_process = [f for f in os.listdir(GENBANK_DIR) if f.endswith((".gb", ".gbk"))]
    print(f"Found {len(files_to_process)} GenBank files to process...")

    with open(OUTPUT_FILE, "w") as out_f:
        out_f.write("SampleID\tTotal_Length_bp\tLSC_Length\tSSC_Length\tIR_Length\tStatus\n")
        
        for filename in sorted(files_to_process):
            filepath = os.path.join(GENBANK_DIR, filename)
            sample_id = os.path.splitext(filename)[0]
            
            total, lsc, ssc, ir, status = get_plastome_structure(filepath)
            
            print(f"Processed: {sample_id} -> Status: {status}")
            
            out_f.write(f"{sample_id}\t{total}\t{lsc}\t{ssc}\t{ir}\t{status}\n")

    print(f"\nProcessing complete! Results saved to {OUTPUT_FILE}")
