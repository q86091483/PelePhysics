#!/usr/bin/env python3
import re
import sys

def find_matching_paren(s, open_idx):
    depth = 0
    i = open_idx
    while i < len(s):
        if s[i] == '(':
            depth += 1
        elif s[i] == ')':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    return -1

def protect_log_calls(text):
    out = []
    i = 0
    n = len(text)
    pattern = re.compile(r'\blog\(')
    while i < n:
        m = pattern.search(text, i)
        if not m:
            out.append(text[i:])
            break
        start = m.start()
        open_paren = m.end() - 1
        close_paren = find_matching_paren(text, open_paren)
        if close_paren == -1:
            raise RuntimeError(f"括号没配对，位置 {start}")

        inner_expr = text[open_paren+1:close_paren]

        if inner_expr.strip().startswith("amrex::max("):
            out.append(text[i:close_paren+1])
            i = close_paren + 1
            continue

        out.append(text[i:start])
        out.append("log(amrex::max(")
        out.append(inner_expr)
        out.append(", 1.e-200))")
        i = close_paren + 1
    return "".join(out)

if __name__ == "__main__":
    infile = sys.argv[1]
    outfile = sys.argv[2]
    start_line = int(sys.argv[3])   # 比如 26235
    end_line = int(sys.argv[4])     # 比如 35891

    with open(infile, "r") as f:
        lines = f.readlines()

    # 行号是1-based，切片要注意+0/-1
    before = lines[:start_line-1]
    target = lines[start_line-1:end_line]
    after = lines[end_line:]

    target_text = "".join(target)
    new_target_text = protect_log_calls(target_text)

    with open(outfile, "w") as f:
        f.writelines(before)
        f.write(new_target_text)
        f.writelines(after)

    print(f"完成，写到 {outfile}")
