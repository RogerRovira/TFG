# MAINTAIN prompt — end-of-session CONTEXT.md update

**What this is:** the recurring prompt for keeping `CONTEXT.md` in sync with the active-nematic solver. Run it at the **end of every working session**, before the context window closes.

**How to use:** attach the **current `.ipynb`** and the **current `CONTEXT.md`**, then paste the prompt block below. (Don't paste the notebook as text — attach the file; pasting has dropped operators before.)

**Companion:** the one-time `INITIALIZE` prompt mints `CONTEXT.md` from scratch; this prompt only updates an existing one.

---

## Prompt

```
Update the attached CONTEXT.md to reconcile it with the current state of my Julia
active-nematic solver before this context window closes. Both the current .ipynb
and the current CONTEXT.md are attached — those are the only ground truth. Ignore
any earlier chat analysis; it may describe a version of the code that no longer
exists.

Section model: the file has DURABLE sections (Model / Numerics & invariants /
Intentional simplifications) and LIVING sections (Current state / Known issues /
Next steps / Data & restart workflow), plus an append-only Session Log.

UPDATE RULES
- DURABLE sections: leave untouched UNLESS the current code now contradicts them.
  If it does, correct them and note the correction in the Session Log. Do not
  re-derive or reword them otherwise.
- LIVING sections: rewrite in place to match the current code. Move completed
  TODOs and now-obvious issues OUT (they survive only as a Session Log line).
- Session Log: append exactly ONE dated entry, newest first — what changed in the
  code, what was decided, what moved on/off the TODO list.
- Bump the "Last updated" line and the code-state pointer.

DISCIPLINE
- Keep it high-signal: add nothing recoverable by reading the code; prune
  ruthlessly. If the file exceeds ~400 lines, compress the oldest Session Log
  entries to one line each.
- Separately, output any NEW inline comments worth adding this session — only on
  non-obvious lines, each tied to a specific line, the why not the what.
```
