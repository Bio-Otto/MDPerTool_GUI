
my_dict = {'anahtar1': 'değer1', 'anahtar2': 'değer2', 'anahtar3': 'değer3', 'anahtar4': 'değer4'}

eski_deger = my_dict.pop('anahtar1')

print(my_dict)



my_dict = {**{'anahtar1': 'yeni_değer1'}, **my_dict}

print(my_dict)

