import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_ga_stats(file_path):
    # Load CSV
    df = pd.read_csv(file_path)

    # Basic sanity check
    if "generation" not in df.columns:
        raise ValueError("The CSV file does not contain a 'generation' column.")

    # Plot Precision and Recall
    plt.figure(figsize=(10, 6))
    plt.plot(df["generation"], df["best_one_precision"], marker="o", label="Precision")
    plt.plot(df["generation"], df["best_one_recall"], marker="s", label="Recall")

    # Optionally plot Avg and Max score
    if "avg" in df.columns and "max" in df.columns:
        plt.plot(df["generation"], df["avg"], linestyle="--", color="gray", label="Average Score")
        plt.plot(df["generation"], df["max"], linestyle=":", color="black", label="Max Score")

    # Title and labels
    plt.title(f"GA Performance Over Generations\n({os.path.basename(file_path)})")
    plt.xlabel("Generation")
    plt.ylabel("Score")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    stats_dir = "stats"

    # List available CSV files
    files = [f for f in os.listdir(stats_dir) if f.endswith(".csv")]
    if not files:
        print(f"No CSV files found in '{stats_dir}/'. Run the GA first.")
        exit()

    print("\nAvailable GA stats files:")
    for i, f in enumerate(files, 1):
        print(f"{i}. {f}")

    # Ask which file to plot
    try:
        choice = int(input("\nEnter the number of the file to plot: ")) - 1
        file_path = os.path.join(stats_dir, files[choice])
    except (ValueError, IndexError):
        print("Invalid choice.")
        exit()

    plot_ga_stats(file_path)
