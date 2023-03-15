from Classes.CancerModel import CancerModel

gridsize     = 201 #PDF: 201
width        = gridsize
height       = gridsize
grids_number = 2

a = CancerModel(388, width, height, grids_number, seed = 10)
print(f'The seed for this test is {a._seed}')
for i in range(48000):
    print(f"Step: {i}")
    a.step()




#Begin travel!
#{1121: [<Classes.CancerCell.CancerCell object at 0x0000027B2EE0E890>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE22350>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE22920>, <Classes.CancerCell.CancerCell object at 0x0000027B2EDF8AF0>]}
#Begin travel!
#{1121: [<Classes.CancerCell.CancerCell object at 0x0000027B2EE0E890>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE22350>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE22920>, <Classes.CancerCell.CancerCell object at 0x0000027B2EDF8AF0>], 1206: [<Classes.CancerCell.CancerCell object at 0x0000027B2EE7CAC0>, <Classes.CancerCell.CancerCell object at 0x0000027B2F198CA0>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE5D0C0>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE5D1E0>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE5FAF0>, <Classes.CancerCell.CancerCell object at 0x0000027B2F19BEB0>, <Classes.CancerCell.CancerCell object at 0x0000027B2F199000>, <Classes.CancerCell.CancerCell object at 0x0000027B2F1B05E0>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE7C7F0>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE5DFC0>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE5D9C0>, <Classes.CancerCell.CancerCell object at 
#0x0000027B2EE417B0>, <Classes.CancerCell.CancerCell object at 0x0000027B2EE22770>]}
#Traceback (most recent call last):
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\DebugTest.py", line 12, in <module>
#    a.step()
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 78, in step
#    self.proliferate("epithelial")
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 89, in proliferate
#    self.schedule.add(new_cell)
#  File "C:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\.venv\lib\site-packages\mesa\time.py", line 68, in add
#    raise Exception(
#Exception: Agent with unique id 7772 already added to scheduler


# Second try got the same results:


#Step: 1119
#Step: 1120
#Begin travel!
#{1121: [<Classes.CancerCell.CancerCell object at 0x0000027593C1A8F0>, <Classes.CancerCell.CancerCell object at 0x0000027593C2E3B0>, <Classes.CancerCell.CancerCell object at 0x0000027593C2E980>, <Classes.CancerCell.CancerCell object at 0x0000027593C0CB50>]}
#Step: 1121
#Step: 1122
#
#Step: 1204
#Step: 1205
#Begin travel!
#{1121: [<Classes.CancerCell.CancerCell object at 0x0000027593C1A8F0>, <Classes.CancerCell.CancerCell object at 0x0000027593C2E3B0>, <Classes.CancerCell.CancerCell object at 0x0000027593C2E980>, <Classes.CancerCell.CancerCell object at 0x0000027593C0CB50>], 1206: [<Classes.CancerCell.CancerCell object at 0x0000027593C90B20>, <Classes.CancerCell.CancerCell object at 0x0000027593FA8D00>, <Classes.CancerCell.CancerCell object at 0x0000027593C69120>, <Classes.CancerCell.CancerCell object at 0x0000027593C69240>, <Classes.CancerCell.CancerCell object at 0x0000027593C6BB50>, <Classes.CancerCell.CancerCell object at 0x0000027593FABF10>, <Classes.CancerCell.CancerCell object at 0x0000027593FA9060>, <Classes.CancerCell.CancerCell object at 0x0000027593FC0640>, <Classes.CancerCell.CancerCell object at 0x0000027593C90850>, <Classes.CancerCell.CancerCell object at 0x0000027593C6A020>, <Classes.CancerCell.CancerCell object at 0x0000027593C69A20>, <Classes.CancerCell.CancerCell object at 
#0x0000027593C55810>, <Classes.CancerCell.CancerCell object at 0x0000027593C2E7D0>]}
#Step: 1206
#Step: 1207


#Step: 1219
#Traceback (most recent call last):
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\DebugTest.py", line 12, in <module>
#    a.step()
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 78, in step
#    self.proliferate("epithelial")
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 89, in proliferate
#    self.schedule.add(new_cell)
#  File "C:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\.venv\lib\site-packages\mesa\time.py", line 68, in add
#    raise Exception(
#Exception: Agent with unique id 7772 already added to scheduler



#Step: 1090
#Step: 1091
#Traceback (most recent call last):
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\DebugTest.py", line 12, in 
#<module>
#    a.step()
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerModel.py", line 83, in step
#    self.schedule.step()
#  File "C:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\.venv\lib\site-packages\mesa\time.py", line 129, in step
#    agent.step()
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerCell.py", line 17, in step
#    self.move()
#  File "c:\Users\vinis\Desktop\Pesquisa IST 2022 - modelagem python células cancer\Repositório-Git\agent-based-cancer\Classes\CancerCell.py", line 51, in move
#    new_position = self.random.choices(possible_steps,weights,k=1)[0]
#  File "C:\Users\vinis\AppData\Local\Programs\Python\Python310\lib\random.py", line 535, in choices
#    raise ValueError('Total of weights must be greater than zero')
#ValueError: Total of weights must be greater than zero