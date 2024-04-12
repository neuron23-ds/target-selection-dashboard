import subprocess

def run_cmd(cmd):
    return subprocess.run(cmd, shell=True, check=True, text=True, capture_output=True)

def mount_bucket(bucket_name, bucket_mount_path):
    cmd = f"mkdir -p {bucket_mount_path}"
    run_cmd(cmd)
    cmd = f"/Users/burkelawlor/go/bin/gcsfuse --implicit-dirs {bucket_name} {bucket_mount_path}"
    run_cmd(cmd)
    print(f"Bucket {bucket_name} mounted at {bucket_mount_path}")

def unmount_bucket(bucket_mount_path):
    cmd = f"umount {bucket_mount_path}"
    run_cmd(cmd)
    cmd = f"rm -r {bucket_mount_path}"
    run_cmd(cmd)
    print(f"Bucket at {bucket_mount_path} unmounted")

