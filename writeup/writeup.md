# Introduction

## Introduction

Frontier models from Anthropic, OpenAI, and Google write quilaity working code for data analysis tasks, but they cost cents to dollars per call resulting a quickly balooning bills. Asking Opus to write a hundred-line bash script every time a postdoc aligns a few FASTQ files spends the budget that should go to harder problems. We ask a simple question: can a frontier model write the recipe once, and a free, small open-weight model running on the lab's own hardware turn that recipe into a working script every time after?

The open-weight model landscape in early 2026 comes from about ten labs:

| Lab | Model family |
|---|---|
| Meta | Llama |
| Alibaba | Qwen |
| DeepSeek | DeepSeek-V4 |
| Mistral | Mistral, Mixtral, Codestral |
| Google | Gemma |
| IBM | Granite |
| Allen AI | OLMo |
| Microsoft | Phi |
| NVIDIA | Nemotron *(new in late 2025)* |
| Cohere | Command, Aya |

Three years ago this list had three names; the open-weight scene now rivals the closed APIs in coverage. Two absences: Apple's on-device 3-billion-parameter model that ships inside iOS and macOS 26 is not downloadable separately, and xAI has released only Grok-1 (March 2024) and Grok-2.5 (August 2025) — every newer Grok is closed.

Sizes span four orders of magnitude. At the small end, Microsoft's Phi-4-mini is 3.8 billion parameters; the 7- to 14-billion class (Phi-4, OLMo 3 7B, Granite 8B, Mistral Ministral) runs comfortably on a single consumer GPU; the 27- to 35-billion dense tier (Qwen3.6-27B, Gemma 4 31B, OLMo 3 32B, Granite 30B) needs around 20 GB of VRAM at 4-bit quantization. The largest flagships are 100-billion- to 1.6-trillion-parameter Mixture-of-Experts models (Llama 4 Scout at 109B, Maverick at 400B, DeepSeek-V4-Pro at 1.6T) that require server-class hardware, though only 13 to 49 billion parameters are active per token. Every flagship released since late 2025 is a sparse Mixture-of-Experts; below 40 billion, dense models persist because they are simpler to deploy, and Alibaba's April 2026 release of a 27-billion dense Qwen3.6 that beats their own 397-billion MoE on coding benchmarks suggests that training quality, not parameter count, now sets the ceiling. Most families also publish a coder-tuned variant (DeepSeek-Coder, Qwen3-Coder-Next, Codestral, Granite-Code) and a "reasoning" or "thinking" variant (DeepSeek-V4-Pro's Max mode, OLMo 3-Think, Phi-4-reasoning, the Qwen3 thinking toggle), an idea from DeepSeek-R1 (January 2025) that every major lab has since copied.

Two practical filters narrow the choice for most labs. First, license: most weights now ship under Apache 2.0 or MIT (DeepSeek, Qwen, Gemma 4, Mistral Large 3, Granite, OLMo, Phi), but Meta's Llama 4 retains a custom Community License with a 700-million-monthly-active-user cap and an explicit exclusion of EU users, and Cohere's Command and Aya are research-and-non-commercial only. Anyone publishing a method built on a given model has to know which license applies. Second, hardware: a 7-billion model at Q4 quantization fits in 5–6 GB of VRAM (any consumer GPU); a 27- to 32-billion dense model fits in a 24 GB GPU (RTX 4090, A5000); the 100-billion-class MoE flagships need 50–60 GB just for weights, which on commodity hardware means two GPUs at minimum; the trillion-parameter tier is multi-node regardless of quantization. Allen AI's OLMo and NVIDIA's Nemotron also release the training data and recipes, not just the weights — relevant when independent replication, not just inference, matters.
