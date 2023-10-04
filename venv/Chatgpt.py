import tkinter as tk
from tkinter import ttk, messagebox


class HydrogeologicalCalculator:
    def __init__(self, root):
        self.root = root
        self.root.title("Hydrogeological Calculator")

        # Create input fields
        self.K_label = ttk.Label(root, text="Hydraulic conductivity (K):")
        self.K_entry = ttk.Entry(root)
        self.vdop_label = ttk.Label(root, text="Vertical Dispersivity (vdop):")
        self.vdop_entry = ttk.Entry(root)
        self.p_label = ttk.Label(root, text="Effective Porosity (p):")
        self.p_entry = ttk.Entry(root)

        # Create calculate button
        self.calculate_button = ttk.Button(root, text="Calculate", command=self.calculate)

        # Create output label
        self.result_label = ttk.Label(root, text="")

        # Place widgets using grid layout
        self.K_label.grid(row=0, column=0, padx=10, pady=5, sticky="w")
        self.K_entry.grid(row=0, column=1, padx=10, pady=5)
        self.vdop_label.grid(row=1, column=0, padx=10, pady=5, sticky="w")
        self.vdop_entry.grid(row=1, column=1, padx=10, pady=5)
        self.p_label.grid(row=2, column=0, padx=10, pady=5, sticky="w")
        self.p_entry.grid(row=2, column=1, padx=10, pady=5)
        self.calculate_button.grid(row=3, columnspan=2, padx=10, pady=10)
        self.result_label.grid(row=4, columnspan=2, padx=10, pady=5)

    def calculate(self):
        try:
            # Get input values
            K = float(self.K_entry.get())
            vdop = float(self.vdop_entry.get())
            p = float(self.p_entry.get())

            # Perform calculation
            result = K * vdop * p

            # Display the result
            self.result_label.config(text=f"Result: {result:.2f}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input. Please enter valid numbers.")


def main():
    root = tk.Tk()
    app = HydrogeologicalCalculator(root)
    root.mainloop()


if __name__ == "__main__":
    main()
