import subprocess

def bash(cmd):
  result = subprocess.run(cmd, shell=True, text=True, capture_output=True)
  if result.returncode != 0:
    raise Exception(f"Command `{cmd}` returned non-zero code {result.returncode}: {result.stderr}")
  return result.stdout


