import subprocess
import sys
import os

def run_script(script_path):
    """Runs a python script and halts execution if it fails."""
    print(f"⏳ Starting: {script_path}...")
    try:
      
        subprocess.run([sys.executable, script_path], check=True)
        print("Finished: {script_path}\n")
    except subprocess.CalledProcessError as e:
        print(" error: {script_path} failed with exit code {e.returncode}.")
        print("Pipeline stopped.")
        sys.exit(1)

def main():
    
    pipeline = [
        os.path.join("data", "deal.py"), # 1. Initial data processing
        "get_paired_lnc_mirna.py",       # 2. Pair generation
        "get_positive_negative.py",      # 3. Positive/negative sample selection
        "train_doc2vec.py",              # 4. Train graph/sequence embeddings
        "get_features.py",               # 5. Extract final features
        "train.py"                       # 6. Train the final classifier
    ]

    print("=== Starting LncMirNet Pipeline ===\n")
    
    for script in pipeline:
        if os.path.exists(script):
            run_script(script)
        else:
            print(f"Missing File: Cannot find '{script}' in the current directory.")
            sys.exit(1)

    print("=== Pipeline Execution Complete ===")

if __name__ == "__main__":
    main()