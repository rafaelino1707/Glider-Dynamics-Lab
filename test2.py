with open("s1223_aero_data.txt", "r") as file:
    i = 0
    lines = file.readlines()
    for line in lines:
        if i >= 12:
            list = line.strip().split()
            print(list)
        i+=1