# Independent-VM experiments

The C client runs ProbGen, Verify, Decode, and native polynomial evaluation. Each server VM runs its own C Compute process. A Python 3 controller transports the binary messages over concurrent TCP connections and records the measurements.

## 1. Prepare the VMs

Use the paper's deployment configuration:

| Role | OS | vCPU / processor | RAM | CPU frequency |
|---|---|---|---|---|
| Client | Ubuntu Desktop 20.04.2 LTS | 1 / Intel i7-6700 | 4 GB | 800 MHz |
| Each server | Ubuntu Desktop 20.04.2 LTS | 1 / Intel Xeon Gold 5218 | 32 GB | 2.30 GHz |

Configure bridged networking and a different hostname and IP address for every VM. The largest server set in the supplied profiles has 11 members. Set CPU allocation and frequency in the virtualization/host configuration before measuring. The runner records the guest's reported CPU, frequency controls, affinity, memory, compiler, libraries, source hash, and every server's metadata.

Copy the same source archive to the client and all servers. Run the following commands from its `Publicly-Verifiable-Outsourcing-of-Polynomial-Computations` directory on every VM:

```sh
sudo apt-get update
sudo apt-get install build-essential libgmp-dev libmpfr-dev m4 curl python3 python3-matplotlib pkg-config
bash scripts/install_paper_dependencies.sh
source .deps/paper/env.sh
make check experiments
```

The installer builds FLINT 2.8.0 and libsodium 1.0.18 under `.deps/paper/`. It verifies the downloaded archives against the pinned SHA-256 hashes. GMP comes from the VM's packages and its loaded version is recorded with each run.

## 2. Start one server on each server VM

On server 0, substitute its bridge IP for the example address:

```sh
source .deps/paper/env.sh
python3 scripts/vm_server.py --server-id 0 --bind 192.168.56.101 --port 9100 > server-0.log 2>&1
```

On subsequent VMs use consecutive IDs 1–10 and their own addresses. Each daemon accepts one active computation at a time. Allow the client VM to connect to TCP port 9100 on the private benchmark network. The transport uses plain TCP on that network.

## 3. Configure the client and run

```sh
source .deps/paper/env.sh
cp experiments/vms.example.json experiments/vms.json
```

Edit the `host` addresses in `experiments/vms.json`. The ordered entries map directly to server IDs. Then run a correctness and connectivity check:

```sh
python3 scripts/run_vm_experiments.py --config experiments/vms.json --profile smoke --trials 3 --output results/measured/vm-smoke
```

Run the figure grids and the client microbenchmarks:

```sh
python3 scripts/run_vm_experiments.py --config experiments/vms.json --profile paper --trials 30 --output results/measured/vm-paper
python3 scripts/run_microbenchmarks.py --trials 100 --output results/measured/vm-micro
python3 scripts/plot_vm_results.py results/measured/vm-paper --output results/plots/vm-paper
```

The protocol repeat count is explicit (`--trials 30` above). The primitive microbenchmarks use the paper's 100 independent trials. Both commands preserve individual samples and report averages and variation.

| Profile | Grid |
|---|---|
| `figure2` | Three schemes; d=2, t=1, m=200,400,…,2000 |
| `figure3` | FE and OV-Shamir; d=2, t=2,3, m=200,400,…,2000 |
| `figure4` | FE and Pi3; d=2,3, t=1,2, m=200,400,…,2000 |
| `figure5` | Three schemes; d=4,8,10, t=1, m=1,…,15 |
| `paper` | Union of Figures 2–5, reusing overlapping configurations |
| `smoke` | All four schemes; d=2, t=1, m=2 |

Use `--profile figure5` to run one grid or `--schemes Gamma_FE_Sh Gamma_OV_Sh` to select schemes. `--dry-run` writes `plan.json` before connecting. `--max-coefficients` sets the dense-coefficient allocation limit on both the client and the server; its default is 5,000,000. `--timeout` and server `--job-timeout` set the per-job deadline in seconds (default 86,400).

Resume an interrupted run by repeating the original command with `--resume`. Completed jobs are retained; unfinished jobs restart. The runner checks the source and run settings before resuming.

## 4. Read the results

Each job directory contains `raw.csv`, `client.stderr`, `servers.json`, and `status.json`. The run directory contains client metadata, the parameter plan, and completion status. A row is written after verification succeeds and Decode equals the native reference.

| Quantity | CSV field / definition |
|---|---|
| T_p | `probgen_ms` |
| T_v | `verify_ms` |
| T_d | `decode_ms` |
| T_c | `client_ms = probgen_ms + verify_ms + decode_ms` |
| T_s | `server_max_ms = max(server_ms)` |
| Per-server work | `server_ms` array in server-ID order; `server_sum_ms` retained for diagnostics |
| T_n | `native_ms`, timed on the client with the same F and x |
| Network round trip | `tcp_roundtrip_ms`, from query upload through the last response |
| Complete online call | `online_wall_ms`, including Python/C IPC, encoding, TCP, and client computation |
| Setup | `keygen_ms`, `remote_setup_ms` |

The controller uploads all queries concurrently, waits for the input-loaded acknowledgements, and dispatches Compute start tokens concurrently. Every server measures Compute with its own monotonic clock. Elapsed intervals need no clock synchronization between VMs. Encoding, parsing, and network wall times have their own columns.

One warmup request precedes recorded trials by default. Public coefficient/input fixtures use SHA-512 labels shared by all schemes; masks, query points, and Paillier keys/coins are freshly sampled. Native evaluation uses an adaptive repetition count targeting 100 ms, capped at 10,000; `native_inner_runs` records the actual count. `--native-repeats N` selects a fixed count. Protocol cleanup occurs after phase timing; native timing includes the evaluator's output cleanup.

Figure CSVs use the ratio of mean client and mean native time, then log10 for Figures 2, 3, and 5. Figure 4 uses mean FE client time divided by mean Pi3 client time for matching parameters. The exporter checks grid coverage, writes `coverage.json`, and produces PDF and PNG figures. `--allow-partial` labels incomplete plots; `--summary-only` exports the measured configurations without figures.

The microbenchmark records Add_q, Mul_q, Add_q2, Mul_q2, random-base Ristretto255 scalar multiplication, and one ordinary randomized Paillier encryption. Each trial samples fresh operands; the first five operations use an inner loop to amortize the clock overhead. Paillier uses one encryption per trial and records key generation separately. Complete HSS Input includes encrypted-digit work and is measured in ProbGen.

## Local verification

```sh
make check experiments
python3 tests/test_vm.py
python3 tests/test_vm.py --hss
```

The TCP test starts temporary localhost workers, checks all decoded results, verifies the maximum-server statistic, rejects malformed frames and wrong server IDs, and checks resume behavior. VM measurement outputs are created by the deployment commands above.
