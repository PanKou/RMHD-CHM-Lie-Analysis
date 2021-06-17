# I want to put the interpretation code I write in a separate file, 
# so I can call them when I need to in a main analysis file


# A function that will interpret the string as a list:

def intlist(s):
    """ Recieves a list written in a string format
        and transforms it into a list object-

        Args:
        s (str):  string version of a list.
    """
    s = s.strip('[')
    s = s.strip(']')
    l = s.split(', ')
    interpreted = []
    for n in l:
        interpreted.append(int(n))
    return interpreted

def get_terms(equation):
    """A function that will split the equations into terms. 
       The terms are always separated by a + or -.

       Args:
       equation (str): equation written in a string form.
    """
    equation = equation.replace('+','-')
    # Break the equation into it's terms
    term_list = equation.split('-')
    # Now strip the whitespace
    terms = [term.strip(' ') for term in term_list]
    if terms[0] == '':
        terms.pop(0)
    # This gets rid of the empty string that results from the equation starting with a negative sign.
    return terms


def get_signs(equation):
    """For a given equation in a string format will give 
       a list with the sign of each term.

       Args:
       equation (str): equation written in a string form.
    """
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


def find_constants(data):
    """Find the constants in the equations in the data set.

       data (DataFrame): dataframe containing all equations.
    """
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
                                                

def process_term(term, cnst_list):
    """ # A function that splits a single term up 
        into individual parts

        Args:
        term (str): Individual term of the equation in
                    string format
        cnst_list (list): list with the constants of the
                          set of original differential equations
    """
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
    fnc, derivatives = process_derivative(parts[-1])
    return [parts[0], greek_list, derivatives, fnc]

# Function to process the derivatives.

def process_derivative(derivative):
    """ Takes the derivative and variable information from 
        a string

        Args:
        derivative (str): A containing the information of which
                          variable is being derivend and the order
                          of the derivatives
    """
    derivative = derivative.strip('Derivative')
    int_list = derivative[:derivative.find(']')+1]
    fnc = derivative[derivative.find(']')+1:].strip('[]')
    return fnc, intlist(int_list)

def term_to_dict(term):
    """ Takes the list with the information of the
        each term and put it in a dictionary 

        Args:
        term (list): List containing the information of
                     an individual term. 
    """
    term_dict = {"coefficient": term[0],
                 "constants": term[1],
                 "derivatives": term[2],
                 "variable": term[3]}
    return term_dict

def eqn_process(equation, cnst_list):
    """ Takes an equation and split it in a list of dictionaries
        saving all the relevant information for each term

        Args:
        term (list): List containing the information of
                     an individual term. 
    """
    signs = get_signs(equation)
    terms = get_terms(equation)
    signs_terms =  zip(signs, terms)
    list_terms = []
    for sign, term in signs_terms:
        processed_term = process_term(term, cnst_list)
        term_dict = term_to_dict(processed_term)
        term_dict["coefficient"] = term_dict["coefficient"]*sign
        list_terms.append(term_dict)
    return list_terms
        