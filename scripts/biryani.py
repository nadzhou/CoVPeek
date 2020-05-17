
from collections import Counter






comment = "trying my best to learn coding so I can apply for Google"

word_bag = comment.split(" ")

word_bag = Counter(word_bag)

print(f"original sentence: {comment}")
print("Word count:")
for k, v in word_bag.items(): 
    print(k, v)