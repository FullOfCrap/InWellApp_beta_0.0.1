#!/usr/bin/python
#!/usr/bin/env python
#-*-coding: utf-8-*-
#STIX - KALKULATRO HYDROGEOLOGICZNY
import math
from tkinter import *
from tkinter import messagebox
import tkinter as tk
from functools import partial
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter.filedialog import askopenfilename
link, warn_list, frame_id, timeStamp = [[] for _ in range(4)]
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import csv
class Welcome():

    def __init__(self, master):
        self.master = master
        self.master.geometry('300x300+100+200')
        self.master.title('Kalkulator STIX - v0.5 Alpha')

        self.label1 = Label(self.master, text='Co chciałbyś zrobić?').grid(row=0)
        self.button1 = Button(self.master, text='Projekt', command=self.projekt, width=40, height=3).grid(row=1)
        self.button2 = Button(self.master, text='Dokumnetacja', command=self.dokumentacja, width=40, height=3).grid(row=2)
        self.button3 = Button(self.master, text='Próbne pompowanie', command=self.pompowanie, width=40, height=3).grid(row=3)

        self.button5 = Button(self.master, text='Wyjdź', command=self.quit, width=40, height=3).grid(row=4)

    def projekt(self):
        root2 = Toplevel(self.master)
        myGUI = projekt(root2)

    def dokumentacja(self):
        root3 = Toplevel(self.master)
        myGUI2 = dokumentacja(root3)

    def pompowanie(self):
        root4 = Toplevel(self.master)
        myGUI3 = pompowanie(root4)

    def quit(self):
        self.master.destroy()


class projekt():

    def __init__(self, master):

        self.wydajnosc = DoubleVar()
        self.srednicaotworu = DoubleVar()
        self.dlugoscfiltra = DoubleVar()
        self.wspolczynnikfiltracj = DoubleVar()
        self.miazszosc = DoubleVar()
        self.parametrliczbowy = DoubleVar()
        self.jednostkowa = DoubleVar()

        self.master = master
        self.master.geometry('700x400+100+200')
        self.master.title('Projekt')

        self.label1 = Label(self.master, text='zakladana wydajnosc studni [m3/h]:').grid(row=0, column=1)
        self.label2 = Label(self.master, text='średnica otworu [m]:').grid(row=1, column=1)
        self.label3 = Label(self.master, text='długość części roboczej filtra [m]:').grid(row=2, column=1)
        self.label4 = Label(self.master, text='współczynnik filtracji w rejonie projektowanego ujęcia [m/h]:').grid(
            row=3, column=1)
        self.label5 = Label(self.master, text='miąższość warstwy wodonośnej [m]: ').grid(row=4, column=1)
        self.label6 = Label(self.master, text='parametr liczbowy dla warstwy wodonosnej w przedziale 1 - 6').grid(row=5, column=1)
        self.label7 = Label(self.master, text='wydajność jednostkowa [m3/h/1mS]:').grid(row=6, column=1)

        self.mywydajnosc = Entry(self.master, textvariable=self.wydajnosc).grid(row=0, column=2)
        self.mysrednicaotworu = Entry(self.master, textvariable=self.srednicaotworu).grid(row=1, column=2)
        self.mydlugoscfiltra = Entry(self.master, textvariable=self.dlugoscfiltra).grid(row=2, column=2)
        self.mywspolcznynnikfiltracji = Entry(self.master, textvariable=self.wspolczynnikfiltracj).grid(row=3, column=2)
        self.mymiazszosc = Entry(self.master, textvariable=self.miazszosc).grid(row=4, column=2)
        self.myparametrliczbowy = Entry(self.master, textvariable=self.parametrliczbowy).grid(row=5, column=2)
        self.myjednostkowa = Entry(self.master, textvariable=self.jednostkowa).grid(row=6, column=2)

        self.button1 = Button(self.master, text='Powierzchnia flitracji [m2]', command=self.powP, width=55).grid(row=7, column=1)
        self.button2 = Button(self.master, text='Dopuszczalna predkosc wlotowa do filtra [m/d]:', command=self.prdop, width=55).grid(row=8, column=1)
        self.button3 = Button(self.master, text='Dopuszczalna wydajność studni [m3/h]:', command=self.qsdop, width=55).grid(row=9, column=1)
        self.button4 = Button(self.master, text='Przewodnosc warstwy [m3/h]:', command=self.Tw, width=55).grid(row=10, column=1)
        self.button5 = Button(self.master, text='Optymalna wydajność studni [m3/h]:', command=self.Qse, width=55).grid(row=11, column=1)
        self.button6 = Button(self.master, text='Projketowana depresja [m]', command=self.Ss, width=55).grid(row=12, column=1)
        self.button7 = Button(self.master, text='Wyjdź', command=self.myquit, width=55).grid(row=13, column=1)

    def powP(self):
        global P
        d = self.srednicaotworu.get()
        l = self.dlugoscfiltra.get()
        P = math.pi * float(l) * float(d)
        self.labelresult=Label(self.master, text='Powierzchnia flitracji: %.2f [m2]' % P).grid(row=7, column=2)

    def prdop(self):
        global vdop
        k = self.wspolczynnikfiltracj.get()
        vdop = (9.8 * math.sqrt(float(k * 24)))
        self.labelresult=Label(self.master, text='dopuszczalna predkosc wlotowa do filtra: %.2f [m/h]' % vdop).grid(row=8, column=2)

    def qsdop(self):
        qsdop = (float(vdop) / 24) * float(P)
        self.labelresult=Label(self.master, text='dopuszczalna wydajność studni: %.2f [m3/h]' % qsdop).grid(row=9, column=2)

    def Tw(self):
        global T
        k = self.wspolczynnikfiltracj.get()
        m = self.miazszosc.get()
        T = k * m
        self.labelresult=Label(self.master, text='przewodnosc warstwy: %2.f [m3/h]' % T).grid(row=10, column=2)

    def Qse(self):
        a = self.parametrliczbowy.get()
        Qe = a * T
        self.labelresult=Label(self.master, text='Wydajność optymalna %.2f [m3/h]' % Qe).grid(row=11, column=2)

    def Ss(self):
        Q = self.wydajnosc.get()
        q = self.jednostkowa.get()
        if q == 0:
            messagebox.showinfo("ERROR", "Nie można dzielić przez 0!!")
        else:
            S = Q/q
            self.labelresult=Label(self.master, text="Pojketowana depresja %.2f [m]" % S).grid(row=12, column=2)


    def myquit(self):
        self.master.destroy()


class dokumentacja():

    def __init__(self, master):

        self.master = master
        self.master.geometry('700x400+100+200')
        self.master.title('Dokumentacja')

        self.number1 = DoubleVar()
        self.number2 = DoubleVar()
        self.number3 = DoubleVar()
        self.number4 = DoubleVar()
        self.number5 = DoubleVar()
        self.number6 = DoubleVar()
        self.number7 = DoubleVar()

        self.labelNum1 = Label(self.master, text='podaj wydajność pompowania pomiarowego (Q) [m3/h]:').grid(row=1, column=0)
        self.labelNum2 = Label(self.master, text='podaj depresję w studni pompowania pomiarowego (Sw) [m]:').grid(row=2, column=0)
        self.labelNum3 = Label(self.master, text='podaj średnicę odwierconego otworu (r) [m]:').grid(row=3, column=0)
        self.labelNum4 = Label(self.master, text='podaj wysykość statycznego zwierciadła wody (H) [m]:').grid(row=4, column=0)
        self.labelNum6 = Label(self.master, text='podaj miąższość warstwy wodonośnej (m) [m]:').grid(row=6, column=0)
        self.labelNum7 = Label(self.master, text='podaj długość części czynnej filtra (L) [m]:').grid(row=7, column=0)

        self.labelResultb = Label(self.master)
        self.labelResultK = Label(self.master)
        self.labelResultR = Label(self.master)
        self.labelResultvdop = Label(self.master)
        self.labelResultp = Label(self.master)
        self.labelResultqdop = Label(self.master)
        self.labelResultqj = Label(self.master)
        self.labelResultKa = Label(self.master)
        self.labelResultKab = Label(self.master)

        self.labelResultb.grid(row=8, column=2)
        self.labelResultK.grid(row=9, column=2)
        self.labelResultR.grid(row=10, column=2)
        self.labelResultKa.grid(row=12, column=3)
        self.labelResultvdop.grid(row=12, column=2)
        self.labelResultp.grid(row=13, column=2)
        self.labelResultqdop.grid(row=14, column=2)
        self.labelResultqj.grid(row=15, column=2)
        self.labelResultKab.grid(row=14, column=3)

        self.entryNum1 = Entry(self.master, textvariable=self.number1).grid(row=1, column=2)
        self.entryNum2 = Entry(self.master, textvariable=self.number2).grid(row=2, column=2)
        self.entryNum3 = Entry(self.master, textvariable=self.number3).grid(row=3, column=2)
        self.entryNum4 = Entry(self.master, textvariable=self.number4).grid(row=4, column=2)
        self.entryNum6 = Entry(self.master, textvariable=self.number6).grid(row=6, column=2)
        self.entryNum7 = Entry(self.master, textvariable=self.number7).grid(row=7, column=2)

        self.kzwssz = partial(self.kzwssz, self.labelResultK, self.labelResultR, self.labelResultb, self.number1, self.number2, self.number3, self.number4, self.number6, self.number7)
        self.kzwssn = partial(self.kzwssn, self.labelResultK, self.labelResultR, self.labelResultb, self.number1, self.number2, self.number3, self.number4, self.number6, self.number7)
        self.kzwnsz = partial(self.kzwnsz, self.labelResultK, self.labelResultR, self.labelResultb, self.number1, self.number2, self.number3, self.number4, self.number6, self.number7)
        self.kzwnsn = partial(self.kzwnsn, self.labelResultK, self.labelResultR, self.labelResultb, self.number1, self.number2, self.number3, self.number4, self.number6, self.number7)
        self.vdopf = partial(self.vdopf, self.labelResultKa, self.labelResultvdop)
        self.pf = partial(self.pf, self.labelResultp, self.number3, self.number7)
        self.qdopf = partial(self.qdopf, self.labelResultKab, self.labelResultqdop)
        self.qjf = partial(self.qjf, self.labelResultqj, self.number1, self.number2)
        self.clear = partial(self.clear, self.labelResultK, self.labelResultR, self.labelResultb, self.labelResultvdop, self.labelResultp, self.labelResultqdop, self.labelResultqj, self.labelResultKa, self.labelResultKab)
        self.buttonCal = Button(self.master, text='współczynnik filtracji dla zwierciadła swobodnego w studni zupełnej', command=self.kzwssz, width=55).grid(row=8, column=0)
        self.buttonCal1 = Button(self.master, text='współczynnik filtracji dla zwierciadła swobodnego w studni niezupełnej', command=self.kzwssn, width=55).grid(row=9, column=0)
        self.buttonCal2 = Button(self.master, text='współczynnik filtracji dla zwierciadła napiętego w studni zupełnej', command=self.kzwnsz, width=55).grid(row=10, column=0)
        self.buttonCal3 = Button(self.master, text='współczynnik filtracji dla zwierciadła napiętego w studni niezupełnej', command=self.kzwnsn, width=55).grid(row=11, column=0)
        self.buttonCal4 = Button(self.master, text='dopuszczalna prędkość wlotowa wody do filtra (Vdop)', command=self.vdopf, width=55).grid(row=12, column=0)
        self.buttonCal5 = Button(self.master, text='powierzchnia części roboczej filtra (P)', command=self.pf, width=55).grid(row=13, column=0)
        self.buttonCal6 = Button(self.master, text='wydajność dopuszczalna filtra (Qdop)', command=self.qdopf, width=55).grid(row=14, column=0)
        self.buttonCal7 = Button(self.master, text='wydajność jednostkowa (q)', command=self.qjf, width=55).grid(row=15, column=0)
        self.buttonCal9 = Button(self.master, text='wyczyść wyniki', command=self.clear, width=55).grid(row=16, column=0)
        self.buttonCal8 = Button(self.master, text='wyjdź', command=self.myquit, width=55).grid(row=20, column=0)

    def kzwssz(self, label_resultK, label_resultR, label_resultb, n1, n2, n3, n4, n6, n7, n8=100):
        global K
        num1 = (n1.get())
        num2 = (n2.get())
        num3 = (n3.get())
        num4 = (n4.get())
        num6 = (n6.get())
        num7 = (n7.get())
        num5 = float(num4) - float(num2)
        num8 = n8
        num9 = float(num3) / 2
        if num5 == 0:
            messagebox.showinfo("ERROR", "Wartość obniżonego zwierciadła wody nie może równać się 0!!")
        elif num5 < 0:
            messagebox.showinfo("ERROR", "Wartość obniżonego zwierciadła wody nie może być mniejsza od 0!!")
        elif float(num3) > float(13):
            messagebox.showinfo("ERROR",
                                "Wartość średnicy odwierconego otworu wynosi %.2f m, czy to aby nie przesada?" % float(
                                    num3))
        else:
            K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            resultK = K
            resultR = R
            label_resultK.config(text='K = %f [m/s]' % resultK)
            label_resultR.config(text='R = %f [m]' % resultR)
            label_resultb.config(text='')

    def kzwssn(self, label_resultK, label_resultR, label_resultb, n1, n2, n3, n4, n6, n7, n8=100):
        global K
        num1 = (n1.get())
        num2 = (n2.get())
        num3 = (n3.get())
        num4 = (n4.get())
        num6 = (n6.get())
        num7 = (n7.get())
        num8 = n8
        num5 = float(num4) - float(num2)
        num9 = float(num3) / 2
        if num5 == 0:
            messagebox.showinfo("ERROR", "Wartość obniżonego zwierciadła wody nie może równać się 0!!")
        elif num5 < 0:
            messagebox.showinfo("ERROR", "Wartość obniżonego zwierciadła wody nie może być mniejsza od 0!!")
        elif float(num5) < float(num7):
            messagebox.showinfo("ERROR",
                                "Wartość obniżonego zwierciadła wody nie może być mniejsza od części czynnej filtra (L)!!")
        elif float(num3) > float(13):
            messagebox.showinfo("ERROR",
                                "Wartość średnicy odwierconego otworu wynosi %.2f m, czy to aby nie przesada?" % float(
                                    num3))
        else:
            b = (math.sqrt(float(num7) / float(num6))) * ((((2 * float(num6)) - float(num7)) / float(num6)) ** 0.25)
            K = (float(num1) * (math.log10(float(num8) / float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2)) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * (math.log10(float(num8)/float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2)) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * (math.log10(float(num8)/float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2)) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * (math.log10(float(num8)/float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2)) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * (math.log10(float(num8)/float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2)) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * (math.log10(float(num8)/float(num9)))) / (
                    (1.366 * ((2 * float(num6)) - float(num2))) * float(num2)) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            resultb = b
            resultK = K
            resultR = R
            label_resultK.config(text='K = %f [m/s]' % resultK)
            label_resultR.config(text='R = %f [m]' % resultR)
            label_resultb.config(text='b = %f' % resultb)

    def kzwnsz(self, label_resultK, label_resultR, label_resultb, n1, n2, n3, n4, n6, n7, n8=100):
        global K
        num1 = (n1.get())
        num2 = (n2.get())
        num3 = (n3.get())
        num4 = (n4.get())
        num6 = (n6.get())
        num7 = (n7.get())
        num5 = float(num4) - float(num2)
        num8 = n8
        num9 = float(num3) / 2
        if num5 == 0:
            messagebox.showinfo("ERROR", "Wartość obniżonego zwierciadła wody nie może równać się 0!!")
        elif num5 < 0:
            messagebox.showinfo("ERROR", "Wartość obniżonego zwierciadła wody nie może być mniejsza od 0!!")
        elif float(num3) > float(13):
            messagebox.showinfo("ERROR",
                                "Wartość średnicy odwierconego otworu wynosi %.2f m, czy to aby nie przesada?" % float(
                                    num3))
        else:
            K = (float(num1) * math.log10(float(num8) / float(num9))) / (2.73 * float(num6) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = (float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            resultK = K
            resultR = R
            label_resultK.config(text='K = %f [m/s]' % resultK)
            label_resultR.config(text='R = %f [m]' % resultR)
            label_resultb.config(text='')

    def kzwnsn(self, label_resultK, label_resultR, label_resultb, n1, n2, n3, n4, n6, n7, n8=100):
        global K
        num1 = (n1.get())
        num2 = (n2.get())
        num3 = (n3.get())
        num4 = (n4.get())
        num6 = (n6.get())
        num7 = (n7.get())
        num5 = float(num4) - float(num2)
        num8 = n8
        num9 = float(num3) / 2
        if num5 == 0:
            messagebox.showinfo("ERROR", "Wartość obniżonego zwierciadła wody nie może równać się 0!!")
        elif num5 < 0:
            messagebox.showinfo("ERROR", "Wartość obniżonego zwierciadła wody nie może być mniejsza od 0!!")
        elif float(num5) < float(num7):
            messagebox.showinfo("ERROR",
                                "Wartość obniżonego zwierciadła wody nie może być mniejsza od części czynnej filtra (L)!!")
        elif float(num3) > float(13):
            messagebox.showinfo("ERROR",
                                "Wartość średnicy odwierconego otworu wynosi %.2f m, czy to aby nie przesada?" % float(
                                    num3))
        else:
            b = (math.sqrt(float(num7) / float(num6))) * ((((2 * float(num6)) - float(num7)) / float(num6)) ** 0.25)
            K = ((float(num1) * math.log10(float(num8) / float(num9))) / (2.73 * float(num6) * float(num2))) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = ((float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = ((float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = ((float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = ((float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            print(K)
            print(R)
            K = ((float(num1) * math.log10(float(R) / float(num9))) / (2.73 * float(num6) * float(num2))) * (1 / b)
            K = K / 3600
            R = 3000 * float(num2) * math.sqrt(K)
            resultb = b
            resultK = K
            resultR = R
            label_resultK.config(text='K = %f [m/s]' % resultK)
            label_resultR.config(text='R = %f [m]' % resultR)
            label_resultb.config(text='b = %f' % resultb)

    def vdopf(self, label_resultKa, label_resultvdop):
        global K
        global vdop
        vdop = (19.6 * math.sqrt(float(K * 3600 * 24)))
        print(vdop)
        resultKa = K
        resultvdop = vdop
        label_resultKa.config(text='(K = %f [m/s])' % resultKa, font=('Verdana', 7), fg='blue')
        label_resultvdop.config(text='Vdop = %f [m/d]' % resultvdop)

    def pf(self, label_resultp, n3, n7):
        global p
        num3 = (n3.get())
        num7 = (n7.get())
        p = math.pi * float(num3) * float(num7)
        print(p)
        resultp = p
        label_resultp.config(text='P = %f [m2]' % resultp)

    def qdopf(self, label_resultKab, label_resultqdop):
        global vdop
        global K
        global p
        qdop = (float(vdop) / 24) * float(p)
        print(qdop)
        resultKab = K
        resultqdop = qdop
        label_resultKab.config(text="K = %f [m/s]" % resultKab, font=('Verdana', 7), fg='blue')
        label_resultqdop.config(text='Qdop %f [m3/h]' % resultqdop)

    def qjf(self, label_resultqj, n1, n2):
        num1 = (n1.get())
        num2 = (n2.get())
        qj = float(num1) / float(num2)
        print(qj)
        resultqj = qj
        label_resultqj.config(text='q = %f [m3/h/1mS]' % resultqj)

    def clear(self, label_resultK, label_resultR, label_resultb, label_resultvdop, label_resultp, label_resultqdop,
              label_resultqj, label_resultKa, label_resultKab):
        global K
        K = 0
        global R
        R = 0
        global vdop
        vdop = 0
        global p
        p = 0
        global qdop
        qdop = 0
        global qj
        qj = 0
        label_resultKa.config(text='')
        label_resultK.config(text='')
        label_resultR.config(text='')
        label_resultb.config(text='')
        label_resultvdop.config(text='')
        label_resultp.config(text='')
        label_resultqdop.config(text='')
        label_resultqj.config(text='')
        label_resultKab.config(text='')



    def myquit(self):
        self.master.destroy()




class pompowanie():

    def __init__(self, master):

        self.master = master
        self.master.geometry('1000x600+100+200')
        self.master.title('Próbne pompowanie')


        menu = Menu(self.master)
        self.master.config(menu=menu)

        file = Menu(menu)

        file.add_command(label='Import', command=self.loadfile)

        file.add_command(label='Wyjdź', command=self.myquit)

        menu.add_cascade(label='Select an option:', menu=file)

        label = Label(self.master, text="Welcome", foreground='purple', font=("Times 20 bold italic"))

        label.pack()

        frame1 = tk.LabelFrame(self.master, labelanchor=NW, height=500, width=500, text='Static Information')

        frame1.pack(fill=BOTH, expand=True)

        text_static = Text(frame1, width=45, height=15, bg='lightgray')

        text_static.pack(side=LEFT, fill=BOTH)

        plot = self.graph()

        self.canvas = FigureCanvasTkAgg(plot, frame1)

        self.canvas.get_tk_widget().pack()

        self.canvas._tkcanvas.pack()

    def loadfile(self):

        global frame_id  # reach outside scope to use frame_id

        filename = askopenfilename(parent=self.master,
                                   filetypes=(("Text File", "*.txt"),
                                              ("All Files", "*.*")),
                                   title='Choose a file')
        with open(filename, 'r')as csvfile:
            plots = csv.reader(csvfile, delimiter=",")
            for row in plots:
                x.append(int(row[1]))
                y.append(int(row[2]))
        newplot = self.graph()
        newplot.canvas = self.canvas
        self.canvas.figure = newplot
        self.canvas.draw()

    def graph(self):
        fig = plt.Figure()
        x = []
        y = []
        ax = fig.add_subplot(111)
        ax.bar(x, y, width=0.5, color='lightgreen')
        return fig

    def myquit(self):
        self.master.destroy()




def main():
    root = Tk()
    myGUIWelcome = Welcome(root)
    root.mainloop()


if __name__ == '__main__':
    main()

