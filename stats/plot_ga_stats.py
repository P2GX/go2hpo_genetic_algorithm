import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_ga_stats(file_path):
    # Load CSV, skipping metadata comments (lines starting with #)
    df = pd.read_csv(file_path, comment="#")

    if "generation" not in df.columns:
        raise ValueError("The CSV file does not contain a 'generation' column.")

    # Create the output folder for plots if missing
    plots_dir = os.path.join(os.path.dirname(file_path), "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # Build plot filename (same name as CSV, but .png)
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    plot_path = os.path.join(plots_dir, f"{base_name}.png")

    # Plot Precision, Recall, Avg, Max
    plt.figure(figsize=(10, 6))
    plt.plot(df["generation"], df["best_one_precision"], marker="o", label="Precision")
    plt.plot(df["generation"], df["best_one_recall"], marker="s", label="Recall")

    if "avg" in df.columns and "max" in df.columns:
        plt.plot(df["generation"], df["avg"], linestyle="--", color="gray", label="Average Score")
        plt.plot(df["generation"], df["max"], linestyle=":", color="black", label="Max Score")

    plt.title(f"GA Performance Over Generations\n({os.path.basename(file_path)})")
    plt.xlabel("Generation")
    plt.ylabel("Score")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    # Save and show
    plt.savefig(plot_path, dpi=300)
    print(f"✅ Plot saved to: {plot_path}")
    plt.show()


if __name__ == "__main__":
    stats_dir = "stats"

    # List all CSV files in stats/
    files = [f for f in os.listdir(stats_dir) if f.endswith(".csv")]
    if not files:
        print(f"No CSV files found in '{stats_dir}/'. Run the GA first.")
        exit()

    print("\nAvailable GA stats files:")
    for i, f in enumerate(files, 1):
        print(f"{i}. {f}")

    try:
        choice = int(input("\nEnter the number of the file to plot: ")) - 1
        file_path = os.path.join(stats_dir, files[choice])
    except (ValueError, IndexError):
        print("Invalid choice.")
        exit()

    plot_ga_stats(file_path)
