# Introduction

## Introduction

Frontier models from Anthropic, OpenAI, and Google write quilaity working code for data analysis tasks, but they cost cents to dollars per call resulting a quickly balooning bills. Asking Opus to write a hundred-line bash script every time a postdoc aligns a few FASTQ files spends the budget that should go to harder problems. We ask a simple question: can a frontier model write the recipe once, and a free, small open-weight model running on the lab's own hardware turn that recipe into a working script every time after?

The open-weight model landscape in early 2026 comes from about ten labs:

| Lab | Family / latest | Sizes (params; active for MoE) | Architecture | Coder variant | Reasoning variant | License |
|---|---|---|---|---|---|---|
| Meta | Llama 4 | Scout 109B (17B-A); Maverick 400B (17B-A); Behemoth ~2T (unreleased) | Sparse MoE; multimodal | — | — | Llama 4 Community License (700M-MAU cap, no EU users) |
| Alibaba | Qwen3.6 | 0.6B → 397B (17B-A); 27B dense | Dense + MoE; thinking toggle | Qwen3-Coder-Next | Qwen3 thinking mode | Apache 2.0 |
| DeepSeek | DeepSeek-V4 | V4-Flash 284B (13B-A); V4-Pro 1.6T (49B-A) | Sparse MoE | DeepSeek-Coder | V4-Pro Max mode | MIT |
| Mistral | Mistral Large 3, Mixtral, Codestral, Ministral | 3B → 675B (41B-A) | Granular MoE flagship; dense small | Codestral 25.08 | — | Apache 2.0 (Codestral non-production) |
| Google | Gemma 4 | E2B, E4B, 26B-MoE, 31B-dense | Dense + MoE; multimodal | — | — | Apache 2.0 |
| IBM | Granite 4.1 / 4.0 | 350M, 1B, 3B, 8B, 30B; 4.0-H-Small 32B (9B-A) | Dense (4.1); hybrid Mamba-2/MoE (4.0) | Granite-Code | — | Apache 2.0 |
| Allen AI | OLMo 3 | 7B, 32B | Dense | — | OLMo 3-Think | Apache 2.0; **fully open** (data + recipes) |
| Microsoft | Phi-4 | mini 3.8B, multimodal 5.6B, 14B | Dense; multimodal | — | Phi-4-reasoning 14B | MIT |
| NVIDIA | Nemotron-3 *(new in late 2025)* | Nano-Omni 30B (3B-A); Super 120B (12B-A) | Hybrid Mamba-Transformer MoE | — | built-in agentic stack | NVIDIA Open Model License; data + recipes |
| Cohere | Command A, Aya | Aya 3.35B → Command A 111B | Dense | — | — | CC-BY-NC 4.0 — **research only** |

Three years ago this list had three names. Two absences: Apple's on-device 3-billion-parameter model that ships inside iOS and macOS 26 is not downloadable separately, and xAI has released only Grok-1 (March 2024) and Grok-2.5 (August 2025) — every newer Grok is closed. Every flagship since late 2025 is a sparse Mixture-of-Experts, but dense persists below ~40 B because it is simpler to deploy, and Alibaba's April 2026 27 B dense Qwen3.6 beats their own 397 B MoE on coding benchmarks — training quality, not parameter count, now sets the ceiling. The reasoning-variant column reflects an idea from DeepSeek-R1 (January 2025) that every major lab has since copied.

Two practical filters narrow the choice for most labs. First, license: most weights now ship under Apache 2.0 or MIT (DeepSeek, Qwen, Gemma 4, Mistral Large 3, Granite, OLMo, Phi), but Meta's Llama 4 retains a custom Community License with a 700-million-monthly-active-user cap and an explicit exclusion of EU users, and Cohere's Command and Aya are research-and-non-commercial only. Anyone publishing a method built on a given model has to know which license applies. Second, hardware: a 7-billion model at Q4 quantization fits in 5–6 GB of VRAM (any consumer GPU); a 27- to 32-billion dense model fits in a 24 GB GPU (RTX 4090, A5000); the 100-billion-class MoE flagships need 50–60 GB just for weights, which on commodity hardware means two GPUs at minimum; the trillion-parameter tier is multi-node regardless of quantization. Allen AI's OLMo and NVIDIA's Nemotron also release the training data and recipes, not just the weights — relevant when independent replication, not just inference, matters.
