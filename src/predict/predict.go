package predict

import (
    "os"
    "os/exec"
    "fmt"
)

func Predict(tmp_cov_output string, output_file string, python_path string) {
	fmt.Println("Predict species...")
    cmd := exec.Command(python_path, "predict/predict.py", tmp_cov_output, output_file)
    cmd.Stdout = os.Stdout
    cmd.Stderr = os.Stderr
    cmd.Run()
    fmt.Println("Saving prediction to:", output_file)
    cmd = exec.Command("rm", tmp_cov_output)
    cmd.Run()
}