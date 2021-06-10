# I want to put the interpretation code I write in a separate file, so I can call them when I need to in a main analysis file

import numpy as np
import pandas as pd

# A function that will interpret the string as a list:

def intlist(s):
    s = s.strip('[')
    s = s.strip(']')
    l = s.split(', ')
    interpreted = []
    for n in l:
        interpreted.append(int(n))
    return interpreted

# A function that will split the equations into terms. The terms are always separated by a + or -.

def get_terms(equation):
    equation = equation.replace('+','-')
    # Break the equation into it's terms
    term_list = equation.split('-')
    # Now strip the whitespace
    terms = [term.strip(' ') for term in term_list]
    if terms[0] == '':
        terms.pop(0)
    # This gets rid of the empty string that results from the equation starting with a negative sign.
    return terms

# A function to get the signs

def get_signs(equation):
    sign_array = []
    if equation[0] != "-":
        sign_array.append(1)
    for char in equation:
        if char == "+":
            sign_array.append(1)
        elif char == '-':
            sign_array.append(-1)
        else:
            continue
    return sign_array

# A function that finds all of the constants

def find_constants(data):
    cnst_list = []
    for equation in data:
        for term in get_terms(equation):
            parts = term.split('*')
            if len(parts) > 1:
                try:
                    parts[0] = int(parts[0])
                except:
                # Here we catch any greek letters, and then put an integer of 1 at the beginning of the list
                    parts[0] = parts[0].strip('(')
                    parts = [1] + parts
            
            if len(parts) > 2:
                # At this point we definitely have a greek letter: print(parts[1:-1])
                for letter in parts[1:-1]:
                    if '^' in letter:
                        # If the constant is raised to a power
                        nc = letter.split('^')[0]
                    else:
                        nc = letter
                    if nc in cnst_list:
                        continue
                    else:
                        cnst_list.append(nc)
    return cnst_list
                                                

# A function that splits a single term up into individual parts

def process_term(term, cnst_list):
    # First need to separate coefficients
    parts = term.split('*')
    # Make the numerical coefficient into an integer
    if len(parts) > 1:
        try:
            parts[0] = int(parts[0])
        except:
        # Here we catch any greek letters, and then put an integer of 1 at the beginning of the list
            parts[0] = parts[0].strip('(')
            parts = [1] + parts

    # Need to strip in 2 stages as to not delete the function information
    parts[-1] = parts[-1].strip('z1, z2, z3, z4, z5, z6, z7, z8, z9])').strip('[')
    greek_list = []
    if len(parts) > 2:
        # At this point we definitely have a greek letter: print(parts[1:-1])
        for cnst in cnst_list:
            i = 0
            # In general we have something like 'η^2'. We need to add the power of each greek letter to the list.
            for letter in parts[1:-1]:
                if cnst in letter:
                    if '^' in letter:
                        i += int(letter.split('^')[1])
                    else:
                        i += 1
                else:
                    continue
            greek_list.append(i)          
    else:
        greek_list = [0 for cnst in cnst_list]
        pass
    return [parts[0], greek_list, parts[-1]]

# A more convinient form for analyzing the term lists
process_terms = lambda terms: [process_term(term) for term in terms]

# Function to process the derivatives.

def process_derivative(derivative):
    derivative = derivative.strip('Derivative')
    int_list = derivative[:derivative.find(']')+1]
    fnc = derivative[derivative.find(']')+1:].strip('[]')
    return fnc, intlist(int_list)