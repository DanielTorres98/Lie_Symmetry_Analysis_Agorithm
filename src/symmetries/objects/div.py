from copy import deepcopy

class Div():
    
    def __init__(self, numerator, denominator) -> None:
        self.numerator = numerator
        self.denominator = denominator

    def __eq__(self, other):
        if isinstance(other, Div):
            num = self.numerator==other.numerator
            den = self.denominator==other.denominator

            return (num and den)
        
        else: 
            return False

    def __repr__(self):
        try: 
            assert len(self.denominator.terms)
            num = f'({self.denominator})'
        except:
            num = f'{self.denominator}'

        return f'{self.numerator}/{num}'
