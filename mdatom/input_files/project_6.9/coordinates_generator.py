n = int(input())

arr = [0, 0, 0]

file = open("coords.inp", "w")
file.write("TITLE")
file.write(n)
for x in range(n):
     for i in range(3):
         file.write(str(arr[i]))
         file.write('   ')
     arr[0] += 0.08814
     if x % 2 == 0:
        arr[2] += 0.1245
     else: 
        arr[2]  = arr[2] - 0.1245
     file.write("\n")   

file.close()
