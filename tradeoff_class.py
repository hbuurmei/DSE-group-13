import numpy as np
import copy
import multiprocessing as mp
import platform

class color:
	def __init__(self, code, name):
		self.HTML = code
		self.name = name
	
class param:
	def __init__(self, name, weight, var, func="LRTS", direc="HB", p=1, Limitype ="minmax", Limit_val=2, roundoff=3):
		self.name = name
		self.func = func
		self.dir = direc
		self.p = p
		self.weight = weight
		self.Ltype = Limitype
		self.val_in = []
		self.val_out = []
		self.l_val = Limit_val
		self.r = roundoff
		self.var = var
	
	def stat(self):
		self.sd = np.std(self.val_in)
		self.mu = np.average(self.val_in)

		if self.Ltype == "minmax":
			self.Lv, self.Hv = min(self.val_in), max(self.val_in)
		elif self.Ltype == "SD":
			if type(self.l_val) != int and type(self.l_val) != float:
				raise Exception("for SD Limits Limit_val must be of type float or int")
			self.Lv, self.Hv = self.mu-self.l_val*self.sd, self.mu+self.l_val*self.sd
		elif self.Ltype == "fixed":
			if len(self.l_val) != 2:
				raise Exception("for fixed Limits Limit_val must contain two values")
			self.Lv, self.Hv = self.l_val[0], self.l_val[1]
		else:
			raise Exception("not valid boundary determination method")

	def func_eval(self, evalv):
		if (evalv <= self.Lv and self.dir == "HB") or (evalv >= self.Hv and self.dir == "LB") :
			return 0
		if (evalv >= self.Hv and self.dir == "HB") or (evalv <= self.Lv and self.dir == "LB") :
			return 1

		if self.dir == "LB":
			temHv, temLv = self.Lv, self.Hv
		else:
			temHv, temLv = self.Hv, self.Lv

		if self.func == "LRTS":
			return (evalv - temLv) / (temHv - temLv)
		elif self.func == "IRTS":
			return (1 - np.exp(-(evalv - temLv) / self.p)) / (1 - np.exp(-(temHv - temLv) / self.p))
		elif self.func == "DRTS":
			return (1 - np.exp(-(temHv - evalv) / self.p)) / (1 - np.exp(-(temHv- temLv) / self.p))
		else:
			raise Exception("not valid scoring scheme")
	
	def set_colors(self, color_list):
		self.color = []
		for val in self.val_out:
			#if self.Limitype != "fixed"
			if val != 1:
				self.color.append(color_list[int(val * len(color_list))])
			else:
				self.color.append(color_list[int(val * len(color_list)) - 1])
			#else:
			#	self.color.append(color_list[int(val * len(color_list))])
	
	
class design:
	def __init__(self, name, sourcelist):
		self.name = name
		self.sourcelist = sourcelist


class tradeoff:
	def __init__(self, design_list, param_list):
		self.param_list = param_list
		self.design_list = design_list


	def get_tradeoff(self):
		self.total = np.zeros(len(self.design_list))
		for i in range(len(self.param_list)):
			param = self.param_list[i]
			param.val_in = np.array([design.sourcelist[i] for design in self.design_list])
			param.stat()
			param.val_out = np.array([param.func_eval(val) for val in param.val_in])
			self.total += param.val_out*param.weight
		
	
	def get_output(self, language = "python", color_list=[], width=10,rot="hor",caption=""):
		def val_s(number):
			#print("i", number)
			res = "" if number >= 0 else "-"
			number = np.fabs(number)
			if number < 10 ** (-20):
				res += str(0)
			elif number < 999 and number > 0.01:
				res+= str(round(number, 3))
			else:
				res += "{:.2e}".format(number)
			if res[-2:] == ".0":
				res = res[:-2]
			#print("o", res)
			#input()
			return res

		if language == "python":
			for param in self.param_list:
				print(param.name, ", \t actual value:", end="\t", sep="")
				for val in param.val_in:
					print(val, end=", \t")
				print()
				print(param.name, ", \t scaled value:", end="\t", sep="")
				for val in param.val_out:
					print(round(val, param.r), end=", \t")
				print()
			print("\t final value:", end="\t", sep="")
			for val in self.total:
				print(round(val, 3), end=", \t")
			print()

		if language == "latex":
			if len(color_list)==0:
				raise Exception("color_list is mandatory for Latex output")
			for c in color_list:
				print("\\definecolor{to-" + c.name + "}{HTML}{" + c.HTML + "}")
			for param in self.param_list:
				param.set_colors(color_list)
			print()
			print("\\begin{table}[H]")
			print("\centering")
			print("\caption{" +caption+ "}")
			print("\label{tab:tradeoff-x}")
			if rot == "ver":
				print("\\begin{adjustbox}{width=0.7\paperheight, angle=-90}")
			else:
				print("\\begin{adjustbox}{width=\\textwidth, angle=0}")
			output = "\\begin{tabular}{|c|l|"
			for param in self.param_list:
				w = width * param.weight if param.weight > 0.15 else width * 0.15
				output += "p{" + str(round(w, 3)) + "cm}|"
				output += "p{" + str(round(w, 3)) + "cm}|"
			output +="c|}\hline"
			print(output)
			
			output = "\multicolumn{2}{|c|}{\multirow{-2}{*}{}}"
			for param in self.param_list:
				output += "& \multicolumn{2}{c|}{\multirow{-2}{*}{}}"
			print(str(output) + "& \multirow{-4}{*}{}\\\\")

			output = "\multicolumn{2}{|c|}{\multirow{-2}{*}{\\textbf{Criteria}}}"
			for param in self.param_list:
				output += "& \multicolumn{2}{c|}{\multirow{-2}{*}{\\textbf{"+ param.name + ", " + str(round(param.weight*100, 2)) +"\%}}}"
			print(str(output) + "& \multirow{-4}{*}{} \\\\ \cline{1-2}")

			output = "\multicolumn{2}{|l|}{\multirow{-2}{*}{}}"
			for param in self.param_list:
				output += "& \multicolumn{2}{c|}{("+ val_s(param.Lv) + ", " + val_s(param.Hv) + ")}"
			print(str(output) + "& \multirow{-4}{*}{} \\\\")
				
			output = "\multicolumn{2}{|l|}{\multirow{-2}{*}{\\textbf{Design Concept}}}"
			for param in self.param_list:
				output += "& \multicolumn{2}{c|}{"
				if param.dir == "HB":
					output += " High Best"
				else:
					output += " Low Best"
				output += ", $\sigma="+ val_s(param.sd) + "$}"

			print(str(output) + "& \multirow{-4}{*}{\\textbf{Total}} \\\\ \hline")
			for i in range(len(self.design_list)):
				design = self.design_list[i]
				output = "\multicolumn{2}{|c|}{}"
				end_output = ""
				k = 4

				for param in self.param_list:
					output += "   & \cellcolor{to-" + param.color[i].name + "} & \cellcolor{to-" + param.color[i].name + "} " + str(param.color[i].name) + ""
					end_output += " \cline{" + str(k) + "-" + str(k) + "} "
					k += 2
				print(str(output) + " & \\\\" + str(end_output))

				output = "\multicolumn{2}{|c|}{}"
				for param in self.param_list:
					output += "   & \multicolumn{2}{c|}{\cellcolor{to-" + param.color[i].name + "}}"
				print(str(output) + "& \\\\")

				output = "\multicolumn{2}{|c|}{\multirow{-3}{*}{" + str(design.name) + "}}"
				for param in self.param_list:
					if param.Lv == 0 and param.Hv == 1:
						output += "   &\multicolumn{2}{c|}{\multirow{-2}{*}{\cellcolor{to-" + param.color[i].name + "} " + val_s(param.val_in[i]) + "}}"
					else:
						output += "   &\multicolumn{2}{c|}{\multirow{-2}{*}{\cellcolor{to-" + param.color[i].name + "} " + val_s(param.val_in[i]) + " $\\rightarrow$ " + val_s(param.val_out[i]) + "}}"
				print(str(output) + " & \multirow{-3}{*}{$" + str(round(self.total[i], 3)) + "$} \\\\ \hline")

			print("\end{tabular}")
			print("\end{adjustbox}")
			print("\end{table}")

class sensitivity:
	def __init__(self, tradeoff, samples=10000):
		self.tro = tradeoff
		self.n = samples
		self.to_tech = False
		self.to_p = False
		self.to_weights = False
		self.per = None
		self.weights = None

	def addto_technical(self, variation):
		self.to_tech = True
		self.to_tech_var = variation

	def addto_p(self, variation):
		self.to_p = True
		self.to_p_var = variation

	def addto_weights(self):
		self.to_weights = True

	def sens(self, n):
		tro_temp = copy.deepcopy(self.tro)
		if self.to_p:
			for param in tro_temp.param_list:
				param.p = np.random.normal(param.p, self, self.to_p_var)

		if self.to_weights:
			total = 0
			for param in tro_temp.param_list:
				param.weight = np.random.normal(param.weight, param.var)
				param.weight = max(0, min(param.weight, 1))
				total += param.weight
				
			for param in tro_temp.param_list:
				param.weight /= total

		if self.to_tech:
			for design in tro_temp.design_list:
				for i in range(len(design.sourcelist)):
					design.sourcelist[i] = np.random.normal(design.sourcelist[i], self.tro.param_list[i].sd*self.to_tech_var)
			
		tro_temp.get_tradeoff()
		weights = [w.weight for w in tro_temp.param_list]
		ret = np.zeros(len(tro_temp.design_list))
		ret[np.where(tro_temp.total == np.amax(tro_temp.total))] = 1
		return ret

	def get_sens(self):
		if platform.system() == "Linux":
			pool = mp.Pool(mp.cpu_count())
			self.per = pool.map(self.sens, range(self.n))
			self.per = np.sum(self.per, axis=0)/self.n
		elif platform.system() == "Windows":
			self.per = []
			for i in range(self.n):
				self.per.append(self.sens(i))
			self.per = np.sum(self.per, axis=0)/self.n
		
	
	def get_RMS(self):
		self.RMS = np.zeros(len(self.tro.design_list))
		for param in self.tro.param_list:
			self.RMS += np.multiply(param.val_in-param.mu, param.val_in-param.mu)/(param.sd*param.sd)