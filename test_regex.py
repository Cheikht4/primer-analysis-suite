import regex
max_e = 2
target = "NNNCTGCATNNN"
p3 = "CTG"
p5 = "CAT"
r = f"(?e){p3}(?:{p5}){{e<={max_e}}}"
print(f"Regex: {r}")
print("Match exact:", bool(regex.search(r, "NNNCTGCATNNN")))
print("Match err in p5:", bool(regex.search(r, "NNNCTGGATNNN")))
print("Match err in p3:", bool(regex.search(r, "NNNCTACATNNN")))
