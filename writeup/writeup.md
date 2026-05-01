# Introduction

## Introduction

Frontier models from Anthropic, OpenAI, and Google write quality working code for data analysis tasks, but they cost cents to dollars per call, resulting in quickly ballooning bills. Asking Opus to write a hundred-line bash script every time a postdoc aligns a few FASTQ files spends the budget that should go to harder problems. We ask a simple question: can a frontier model write the recipe once, and a free, small open-weight model running on the lab's own hardware turn that recipe into a working script every time after?

To decide which open models to use we first need to survey the landscape of available models. For 2026 this landscape is summarized in Table 1 below. 

**Table 1.** Major open-weight model families as of May 2026.

| Lab | Family / latest | Sizes (total; active for MoE) | Architecture | Coder variant | Reasoning variant | License |
|---|---|---|---|---|---|---|
| Meta | Llama 4 [1] | Scout 109 B (17 B-A); Maverick 400 B (17 B-A); Behemoth ~2 T (unreleased) | Sparse MoE; multimodal | — | — | Llama 4 Community License (700 M-MAU cap; excludes EU users) |
| Alibaba | Qwen3.6 [2] | 0.6 B → 397 B (17 B-A); 27 B dense | Dense + MoE; thinking toggle | Qwen3-Coder-Next | Qwen3 thinking mode | Apache 2.0 |
| DeepSeek | DeepSeek-V4 [3] | V4-Flash 284 B (13 B-A); V4-Pro 1.6 T (49 B-A) | Sparse MoE | DeepSeek-Coder | V4-Pro Max mode | MIT |
| Mistral | Mistral Large 3, Mixtral, Codestral, Ministral [4] | 3 B → 675 B (41 B-A) | Granular MoE flagship; dense small | Codestral 25.08 | — | Apache 2.0 (Codestral non-production) |
| Google | Gemma 4 [5] | E2B, E4B, 26 B-MoE, 31 B-dense | Dense + MoE; multimodal | — | — | Apache 2.0 |
| IBM | Granite 4.1 / 4.0 [6] | 350 M, 1 B, 3 B, 8 B, 30 B; 4.0-H-Small 32 B (9 B-A) | Dense (4.1); hybrid Mamba-2 / MoE (4.0) | Granite-Code | — | Apache 2.0 |
| Allen AI | OLMo 3 [7] | 7 B, 32 B | Dense | — | OLMo 3-Think | Apache 2.0; **fully open** (weights + data + recipes) |
| Microsoft | Phi-4 [8] | mini 3.8 B, multimodal 5.6 B, 14 B | Dense; multimodal | — | Phi-4-reasoning 14 B | MIT |
| NVIDIA | Nemotron-3 *(new in late 2025)* [9] | Nano-Omni 30 B (3 B-A); Super 120 B (12 B-A) | Hybrid Mamba-Transformer MoE | — | built-in agentic stack | NVIDIA Open Model License; data + recipes |
| Cohere | Command A, Aya [10] | Aya 3.35 B → Command A 111 B | Dense | — | — | CC-BY-NC 4.0 — **research only** |

**Legend.**

- **Architecture.** Modern language models read input one piece at a time (each piece, called a *token*, is roughly a short word or word-fragment) and predict the next piece, by passing each token through a deep stack of trainable mathematical operations. Three architectural styles appear in the table:
    - *Dense.* Every internal weight is used on every token — the original transformer design and what most pre-2024 LLMs were. Simple to deploy, but compute per token scales with the model's full size.
    - *Sparse Mixture-of-Experts (MoE).* The model is split into many parallel "expert" sub-networks plus a small "router" that, for each token, picks one or two experts to handle it; the other experts stay idle. This lets the model be very large in *total* parameters (high stored capacity) while spending the per-token compute of a much smaller model. The cost is memory: every expert still has to be loaded, even the ones that won't fire on a given input. Every flagship released since late 2025 uses this design.
    - *Hybrid Mamba-Transformer.* A newer variant that replaces some standard layers with *state-space* (Mamba) layers. The practical consequence is that compute on a long input grows linearly with input length rather than quadratically — relevant now that flagships accept inputs of 500,000 to 1,000,000 tokens (a small book).
    - *Multimodal.* Accepts images, and sometimes audio, alongside text.
- **Parameters (B = billion, T = trillion).** Parameters are the trained numerical weights inside a model — loosely analogous to synapse strengths. More parameters mean more stored capacity but also more memory and more arithmetic per token. Two numbers matter:
    - *Total* parameters set disk and GPU-memory cost. A 70 B-parameter model occupies roughly 35–40 GB of memory at the standard 4-bit compression used for inference (a 16-bit "full precision" version would need ~140 GB).
    - *Active* parameters (notation `17 B-A` = 17 billion active per token) set the per-token compute cost. For dense models the two numbers are equal. For MoE models the active count is much smaller than the total because most experts don't fire on a given token. A 109 B (17 B-A) MoE runs about as fast per token as a 17 B dense model, but the machine still has to keep all 109 B parameters loaded.
- **Coder variant.** A general model further trained on a large corpus of source code (and tested on standard code-task benchmarks such as SWE-bench and HumanEval). Coder variants beat their general siblings on bash, Python, and refactor tasks but lose ground on chat, math, and translation.
- **Reasoning variant.** A model that, before answering, writes out a step-by-step reasoning trace ("chain of thought") that the user typically does not see; training rewards the model when its final answer is verifiably correct. Pioneered by DeepSeek-R1 (January 2025); every major lab has since copied the recipe. Reasoning variants score higher on math, code, and multi-step problems but take roughly 5×–50× longer per call (and cost 5×–50× more).

Table 1 has notable absences: Apple's 3-billion-parameter on-device model ships only inside iOS and macOS 26, and xAI has openly released only Grok-1 (March 2024) and Grok-2.5 (August 2025) — everything newer is closed. Every flagship since late 2025 is a sparse MoE, but dense persists below ~40 B because it is simpler to deploy, and Alibaba's April 2026 27 B dense Qwen3.6 beats their own 397 B MoE on coding benchmarks — training quality, not parameter count, now sets the ceiling.

Two practical filters narrow the choice for most labs. First, license: most weights now ship under Apache 2.0 or MIT (DeepSeek, Qwen, Gemma 4, Mistral Large 3, Granite, OLMo, Phi), but Meta's Llama 4 retains a custom Community License with a 700-million-monthly-active-user cap and an explicit exclusion of EU users, and Cohere's Command and Aya are research-and-non-commercial only. Anyone publishing a method built on a given model has to know which license applies. Second, hardware: the four size tiers in Table 1 map to four very different machines, from a \$400–\$600 consumer card for the smallest models to a multi-GPU server costing several hundred thousand dollars for the trillion-parameter tier (Table 2). Allen AI's OLMo and NVIDIA's Nemotron also release the training data and recipes, not just the weights — relevant when independent replication, not just inference, matters.

**Table 2.** GPU options and ballpark May 2026 US street prices for running each Table 1 model class locally at 4-bit compression. Where multiple cards are needed, the listed price is the per-card cost.

| Model class (Table 1 example) | GPU memory needed | GPU options (May 2026 USD) | Refs |
|---|---|---|---|
| **7–8 B dense** (Phi-4, Granite 8B, OLMo 3 7B) | ~5–6 GB | NVIDIA RTX 5060 Ti 16 GB (~\$560); RTX 4060 Ti 16 GB (~\$430); RTX 5070 12 GB (~\$635); Apple M4 Mac mini 16 GB unified (~\$600) | [11, 16] |
| **27–32 B dense** (Qwen3.6-27B, Gemma 4 31B, OLMo 3 32B) | ~17–22 GB | NVIDIA RTX 5090 32 GB (~\$2,900–\$3,500, street price 50–75 % over \$1,999 MSRP); RTX 4090 24 GB (~\$1,500–\$2,200, EOL Oct 2024); RTX A5000 24 GB (~\$700–\$1,400 used); Apple Mac Studio M3 Ultra 96 GB unified (~\$4,000) | [11, 12, 16] |
| **100–400 B MoE** (Llama 4 Scout 109B, Maverick 400B) | ~55–250 GB | NVIDIA RTX Pro 6000 Blackwell 96 GB (~\$8,500); 2× RTX 5090 (~\$6,000–\$7,000 total); RTX 6000 Ada 48 GB (~\$6,800); H100 80 GB (~\$25K–\$33K); AMD Instinct MI300X 192 GB (~\$15K–\$20K, OEM-only); Apple Mac Studio M3 Ultra 256 GB unified (~\$9,500) | [11, 12, 13, 14, 16] |
| **Trillion-parameter MoE** (DeepSeek-V4-Pro 1.6T) | ~700+ GB | NVIDIA H200 141 GB (~\$31K–\$40K each); B200 192 GB (~\$35K–\$55K each); 8-GPU DGX H200 server (~\$350K–\$500K); GB200 NVL72 rack (~\$3M+); cloud rental on B200 ~\$2.25–\$16/GPU-hour | [13, 15] |

May 2026 prices are dominated by an ongoing HBM/GDDR7 shortage; consumer NVIDIA 50-series cards and high-RAM Apple Mac Studio configurations sit 30–75 % above launch MSRP. AMD MI300X and the newer MI325X are sold almost exclusively through OEM channels — public per-card numbers reflect bulk pricing or cloud rental rates [13]. The used market is liquid for retired flagships (RTX 4090, RTX A5000) and is often the cheapest path into the 27–32 B tier [11].

Given this (rapidly evolving) landscape we decided to do the following experiment: take a common sequencing data processing workflow and ask open models running on hardware accessible to an average research lab to design and execute the analysis. In doing so we experimented with a range of possibilities ranging from allowing open models to figure out everything by themselves to guiding them using a very detailed plan produced by commercial frontier models. We further complicated these tasks by simulating a variety of errors that may occur during workflow execution.

## Results

### Hardware

We selected five different computers listed in Table 3. It is a combination ranging from an old workstation saved from university salvage and a desktop with a gaming GPU to the latest MacBooks and a purpose-built inexpensive NVIDIA device---NVIDIA Jetson AGX Orin. The Orin is a "RaspberryPi"-like offering from NVIDIA that costs under $2,000 and has a very small footprint making it an idea lab-ready tiny but powerful workstation.

**Table 3.** Test machines used in this study.

| Computer | Manufacturer | Year released | RAM | OS | GPU |
|---|---|---|---|---|---|
| NVIDIA Jetson AGX Orin Developer Kit | NVIDIA | 2022 | 64 GB LPDDR5 (unified with GPU) | Ubuntu 22.04 LTS (NVIDIA JetPack 6) | Integrated Ampere, 2,048 CUDA cores + 64 Tensor cores |
| RTX 5080 desktop | Custom build | 2025 (GPU) | 64 GB DDR5 (system) | Linux (Ubuntu) | NVIDIA RTX 5080, 10,752 CUDA cores, 16 GB GDDR7 |
| MacBook Air M4 (24 GB) | Apple | 2025 | 24 GB LPDDR5X (unified with GPU) | macOS 26 | Apple M4 integrated GPU, 10 cores |
| MacBook Pro M4 Pro (48 GB) | Apple | 2024 | 48 GB LPDDR5X (unified with GPU) | macOS Sequoia 15.6 | Apple M4 Pro integrated GPU, 16 to 20 cores |
| 2× NVIDIA RTX A5000 desktop | Custom build | 2021 (GPUs) | 128 GB DDR4 (system) | Linux (Ubuntu 25.10) | 2× NVIDIA RTX A5000, 8,192 CUDA cores each, 24 GB GDDR6 each (48 GB total) |

System RAM listed for the two Linux desktops reflects the build configuration; for inference workloads the relevant memory is the GPU VRAM (last column). For the Jetson and the MacBook, RAM is unified between CPU and GPU and the model can use up to roughly the listed RAM minus the operating-system reservation.

### Workflow




## References

[1] Meta. Llama 4: a new crop of flagship AI models. *TechCrunch*, April 5, 2025. https://techcrunch.com/2025/04/05/meta-releases-llama-4-a-new-crop-of-flagship-ai-models/

[2] Alibaba (Qwen team). Qwen3.6 family. GitHub. https://github.com/QwenLM/Qwen3.6

[3] DeepSeek-AI. DeepSeek-V4 preview release notes. DeepSeek API documentation, April 24, 2026. https://api-docs.deepseek.com/news/news260424

[4] Mistral AI. Introducing Mistral 3. Mistral AI news, December 2, 2025. https://mistral.ai/news/mistral-3

[5] Google. Introducing Gemma 4. Google blog, April 2, 2026. https://blog.google/innovation-and-ai/technology/developers-tools/gemma-4/

[6] IBM Research. Granite 4.1 AI foundation models. April 30, 2026. https://research.ibm.com/blog/granite-4-1-ai-foundation-models

[7] Allen Institute for AI (Ai2). OLMo 3. November 20, 2025. https://allenai.org/blog/olmo3

[8] Microsoft. Welcome to the new Phi-4 models — Phi-4-mini and Phi-4-multimodal. Microsoft TechCommunity, Educator Developer Blog. https://techcommunity.microsoft.com/blog/educatordeveloperblog/welcome-to-the-new-phi-4-models---microsoft-phi-4-mini--phi-4-multimodal/4386037

[9] NVIDIA. NVIDIA debuts Nemotron 3 family of open models. NVIDIA Newsroom. https://nvidianews.nvidia.com/news/nvidia-debuts-nemotron-3-family-of-open-models

[10] Cohere. Models — Command A, Aya. Cohere docs. https://docs.cohere.com/docs/models

[11] BestValueGPU. Consumer NVIDIA RTX GPU price history and specifications (RTX 4060 Ti, 4090, 5060 Ti, 5070, 5090). https://bestvaluegpu.com/

[12] Thunder Compute. NVIDIA RTX Pro 6000 Blackwell pricing analysis. https://www.thundercompute.com/blog/nvidia-rtx-pro-6000-pricing

[13] Thunder Compute. AMD Instinct MI300X pricing. https://www.thundercompute.com/blog/amd-mi300x-pricing

[14] Jarvis Labs. NVIDIA H100 80 GB pricing guide. https://jarvislabs.ai/blog/h100-price

[15] Northflank. NVIDIA B200 cost analysis and cloud rental rates. https://northflank.com/blog/how-much-does-an-nvidia-b200-gpu-cost

[16] Apple. Mac Studio configurations and pricing. https://www.apple.com/mac-studio/specs/
