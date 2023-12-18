# HF_SVAR_Identification
The repository contains replication codes for implementing heterogenous VAR estimation using high-frequency data as inputs. The methodology is taken from the paper "SVAR Identification with High-Frequency Macroeconomic Data". The goal is to achieve the full model identification of the structural shocks by employing high-frequency macroeconomic data that allows a natural diagonalization of the contemporaneous relations.

## The Model: Heterogeous estimation
$B_0\boldsymbol{y}_{t} = B_d \boldsymbol{y}^{(d)}  + B_w \boldsymbol{y}^{(w)} + B_m \boldsymbol{y}^{(m)} + B_q \boldsymbol{y}^{(q)} + \boldsymbol{e}_t,$

with:

$\boldsymbol{y}^{(d)} = \boldsymbol{y}_{t-1},$

$\boldsymbol{y}^{(w)} = \Sigma_{i=2}^{5} \boldsymbol{y}_{t-i},$

$\boldsymbol{y}^{(m)} = \Sigma_{i=6}^{22} \boldsymbol{y}_{t-i},$

$\boldsymbol{y}^{(q)} = \Sigma_{i=23}^66 \boldsymbol{y}_{t-i}.$


## Natural diagonalization
The contemporaneous correlation between the daily (HF) macroeconomic series disappears as the frequency increases, allowing for full identification of the structural shocks.

![corr_mean_data_abs](https://github.com/l-longo/HF_SVAR_Identification/assets/96943070/2c702932-bf4b-4327-93ab-99e9edb40a8d)

## Getting started
