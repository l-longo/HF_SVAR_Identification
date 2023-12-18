# HF_SVAR_Identification
The repository contains replication codes for implementing heterogenous VAR estimation using high-frequency data as inputs. The methodology is taken from the paper "SVAR Identification with High-Frequency Macroeconomic Data". The goal is to achieve the full model identification of the structural shocks by employing high-frequency macroeconomic data that allows a natural diagonalization of the contemporaneous relations.

# The Model

$B_0\bm{y}_{t} = B_d \textbf{y}^{(d)}  + B_w \bm{y}^{(w)} + B_m \bm{y}^{(m)} + B_q \bm{y}^{(q)} + \bm{e}_t,$

where $\bm{y}^{(d)} = \bm{y}_{t-1}$, $\bm{y}^{(w)} = \sum_{i=2}^5 \bm{y}_{t-i}$, $\bm{y}^{(m)} = \sum_{i=6}^{22} \bm{y}_{t-i}$, $\bm{y}^{(q)} = \sum_{i=23}^{66} \bm{y}_{t-i}$. 


