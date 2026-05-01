# Introduction

## Introduction

Frontier models from Anthropic, OpenAI, and Google write quilaity working code for data analysis tasks, but they cost cents to dollars per call resulting a quickly balooning bills. Asking Opus to write a hundred-line bash script every time a postdoc aligns a few FASTQ files spends the budget that should go to harder problems. We ask a simple question: can a frontier model write the recipe once, and a free, small open-weight model running on the lab's own hardware turn that recipe into a working script every time after?

The open-weight model landscape in early 2026 comes from about ten labs (Table 1).

**Table 1.** Major open-weight model families as of May 2026.

| Lab | Family / latest | Sizes (total; active for MoE) | Architecture | Coder variant | Reasoning variant | License |
|---|---|---|---|---|---|---|
| Meta | Llama 4 | Scout 109 B (17 B-A); Maverick 400 B (17 B-A); Behemoth ~2 T (unreleased) | Sparse MoE; multimodal | — | — | Llama 4 Community License (700 M-MAU cap; excludes EU users) |
| Alibaba | Qwen3.6 | 0.6 B → 397 B (17 B-A); 27 B dense | Dense + MoE; thinking toggle | Qwen3-Coder-Next | Qwen3 thinking mode | Apache 2.0 |
| DeepSeek | DeepSeek-V4 | V4-Flash 284 B (13 B-A); V4-Pro 1.6 T (49 B-A) | Sparse MoE | DeepSeek-Coder | V4-Pro Max mode | MIT |
| Mistral | Mistral Large 3, Mixtral, Codestral, Ministral | 3 B → 675 B (41 B-A) | Granular MoE flagship; dense small | Codestral 25.08 | — | Apache 2.0 (Codestral non-production) |
| Google | Gemma 4 | E2B, E4B, 26 B-MoE, 31 B-dense | Dense + MoE; multimodal | — | — | Apache 2.0 |
| IBM | Granite 4.1 / 4.0 | 350 M, 1 B, 3 B, 8 B, 30 B; 4.0-H-Small 32 B (9 B-A) | Dense (4.1); hybrid Mamba-2 / MoE (4.0) | Granite-Code | — | Apache 2.0 |
| Allen AI | OLMo 3 | 7 B, 32 B | Dense | — | OLMo 3-Think | Apache 2.0; **fully open** (weights + data + recipes) |
| Microsoft | Phi-4 | mini 3.8 B, multimodal 5.6 B, 14 B | Dense; multimodal | — | Phi-4-reasoning 14 B | MIT |
| NVIDIA | Nemotron-3 *(new in late 2025)* | Nano-Omni 30 B (3 B-A); Super 120 B (12 B-A) | Hybrid Mamba-Transformer MoE | — | built-in agentic stack | NVIDIA Open Model License; data + recipes |
| Cohere | Command A, Aya | Aya 3.35 B → Command A 111 B | Dense | — | — | CC-BY-NC 4.0 — **research only** |

**Legend.**

- **Parameters (B = billion, T = trillion).** A model's parameters are its trained weights; more parameters generally mean more capability but also more memory and compute per call. *Total* parameters set disk and VRAM cost. For Mixture-of-Experts models, *active* parameters (notation `17 B-A` = 17 billion active per token) set the per-call compute cost, because only a subset of the experts fire on any given token. A 109 B (17 B-A) MoE runs about as fast per token as a 17 B dense model but still needs all 109 B parameters in memory.
- **Architecture.** *Dense* means every parameter is used on every token — the original transformer design. *Sparse Mixture-of-Experts (MoE)* routes each token through a small subset of "expert" sub-networks, trading memory for compute. *Hybrid Mamba-Transformer* replaces some attention layers with state-space (Mamba) blocks, which scale linearly with context length instead of quadratically — useful at the 512 K- to 1 M-token contexts the new flagships now ship. *Multimodal* models accept images (and sometimes audio) alongside text.
- **Coder variant.** A model post-trained on code-heavy data and benchmarked on code tasks (SWE-bench, HumanEval). Coder variants beat their general siblings on bash, Python, and refactor tasks but lose ground on chat, math, and translation.
- **Reasoning variant.** A model that produces an explicit chain-of-thought before its answer, typically trained with reinforcement learning from verifiable rewards. The recipe was popularised by DeepSeek-R1 (January 2025); every major lab has since copied it. Reasoning variants score higher on math, code, and multi-step problems but spend 5×–50× more output tokens per call.

Three years ago this list had three names. Two absences: Apple's 3-billion-parameter on-device model ships only inside iOS and macOS 26, and xAI has released only Grok-1 (March 2024) and Grok-2.5 (August 2025) — everything newer is closed. Every flagship since late 2025 is a sparse MoE, but dense persists below ~40 B because it is simpler to deploy, and Alibaba's April 2026 27 B dense Qwen3.6 beats their own 397 B MoE on coding benchmarks — training quality, not parameter count, now sets the ceiling.

Two practical filters narrow the choice for most labs. First, license: most weights now ship under Apache 2.0 or MIT (DeepSeek, Qwen, Gemma 4, Mistral Large 3, Granite, OLMo, Phi), but Meta's Llama 4 retains a custom Community License with a 700-million-monthly-active-user cap and an explicit exclusion of EU users, and Cohere's Command and Aya are research-and-non-commercial only. Anyone publishing a method built on a given model has to know which license applies. Second, hardware: a 7-billion model at Q4 quantization fits in 5–6 GB of VRAM (any consumer GPU); a 27- to 32-billion dense model fits in a 24 GB GPU (RTX 4090, A5000); the 100-billion-class MoE flagships need 50–60 GB just for weights, which on commodity hardware means two GPUs at minimum; the trillion-parameter tier is multi-node regardless of quantization. Allen AI's OLMo and NVIDIA's Nemotron also release the training data and recipes, not just the weights — relevant when independent replication, not just inference, matters.
