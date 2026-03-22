# postgwas_pleio

## Overview

This module documents key assumptions and requirements for meta-analysis workflows using **METAL** and **MTAG** within the `postgwas_pleio` pipeline.

---

## Key Notes

### 1. METAL-based Meta-analysis

- **Effective Sample Size (Neff)** must be used.
- During METAL-based meta-analysis:
  - ❌ **Sample prevalence** is **not used**
  - ❌ **Population prevalence** is **not used**

---

### 2. MTAG-based Meta-analysis

- **Effective Sample Size (Neff)** must be used.
- During MTAG-based meta-analysis:
  - ❌ **Sample prevalence** is **not used**
  - ❌ **Population prevalence** is **not used**

---

## References

- **ref1**: https://github.com/JonJala/mtag/issues/239  
- **ref2**: https://github.com/JonJala/mtag/issues/60  

---

## Summary

| Method | Uses Neff | Uses Sample Prevalence | Uses Population Prevalence |
|--------|----------|----------------------|----------------------------|
| METAL  | ✅ Yes   | ❌ No                | ❌ No                      |
| MTAG   | ✅ Yes   | ❌ No                | ❌ No                      |

---

## Notes for Implementation

- Always compute **Neff** using:

  ```text
  Neff = 4 / (1/N_case + 1/N_control)