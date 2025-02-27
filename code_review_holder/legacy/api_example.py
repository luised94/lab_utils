# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 20:15:31 2020

@author: liusm
"""

#Example for using APIs. Uses exchange rate and api. Taken from https://python.gotrained.com/python-json-api-tutorial/
import requests
import json

url = "https://api.exchangeratesapi.io/latest?symbols=USD,GBP"

response = requests.get(url)
data = response.text
parsed = json.loads(data)
date = parsed["date"]

gbp_rate = parsed["rates"]["GBP"]
usd_rate = parsed["rates"]["USD"]
print("On " + date + " EUR equals " + str(gbp_rate) + " GBP")
print("On " + date + " EUR equals " + str(usd_rate) + " USD")




url = "https://api.exchangeratesapi.io/latest?base=USD"

response = requests.get(url)
data = response.text
parsed = json.loads(data)
date = parsed["date"]
print("Date:", date, "\n")

rates = parsed["rates"]

for currency, rate in rates.items():   #use items() method to iterate through keys and values of a dictionary 
    print("USD =",currency, rate)

#changing url parameters
    
bases = ["USD", "EUR", "GBP"]

for base in bases:
    url = "https://api.exchangeratesapi.io/latest?base=" + base

    response = requests.get(url)
    data = response.text
    parsed = json.loads(data)

    rates = parsed["rates"]

    print("--------------- Rates in", base, "---------------")
    for currency, rate in rates.items():
        print(base, "=", currency, rate)


#simple currency converter
base = input("Convert from: ")
to = input("Convert to: ")
amount = float(input("Amount: "))

url = "https://api.exchangeratesapi.io/latest?base=" + base

response = requests.get(url) 
data = response.text 
parsed = json.loads(data) 
rates = parsed["rates"]


for currency, rate in rates.items():
    if currency == to:
        conversion = rate * amount
        print("1", base, "=", currency, rate)
        print(amount, base, "=", currency, conversion)

#with access key
url = "http://data.fixer.io/api/latest?access_key=9d1f065a6a9816f9dd7667b0c9b26bab&symbols=USD,GBP"

response = requests.get(url)
data = response.text
parsed = json.loads(data)
date = parsed["date"]

gbp_rate = parsed["rates"]["GBP"]
usd_rate = parsed["rates"]["USD"]
print("On " + date + " EUR equals " + str(gbp_rate) + " GBP")
print("On " + date + " EUR equals " + str(usd_rate) + " USD")
