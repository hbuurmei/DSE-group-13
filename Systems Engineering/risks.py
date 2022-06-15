import pandas as pd
import os


def get_file(file):
    # Returns the full path, if the file is in the same folder as the main .py program.
    return os.path.join(os.path.dirname(file), file)


def get_folder_file(folder, file):
    # Returns the full path, if the file is not in the same folder as the main .py program.
    # If this does not work, use: return get_file(os.path.join(folder, file))
    return os.path.join(folder, file)


def get_colour(l, c):
    # l = likelihood, c = consequence
    if (c == 'Negligible' and l == 'Very Likely') or (c == 'Minor' and (l == 'Likely' or l == 'Possible' or l == "Unlikely"))\
            or (c == 'Moderate' and (l == 'Unlikely' or l == 'Very Unlikely')):
        return r'\cellcolor{lightgreen} '
    elif (c == 'Negligible' and (l == 'Likely' or l == 'Possible' or l == 'Unlikely' or l == 'Very Unlikely')) or \
            (c == 'Minor' and l == 'Very Unlikely'):
        return r'\cellcolor{forestgreen} '
    elif (c == 'Minor' and l == 'Very Likely') or (c == 'Moderate' and (l == 'Likely' or l == 'Possible')) or \
            (c == 'Significant' and (l == 'Unlikely' or l == 'Very Unlikely')) or (c == 'Severe' and l == 'Very Unlikely'):
        return r'\cellcolor{lightyellow} '
    elif (c == 'Moderate' and l == 'Very Likely') or (c == 'Significant' and (l == 'Likely' or l == 'Possible')) or \
            (c == 'Severe' and (l == 'Possible' or l == 'Unlikely')):
        return r'\cellcolor{lightorange} '
    else:
        return r'\cellcolor{lightred} '


# ----------------------------- Read csv file --------------------------------
total = pd.read_csv(get_file('risks.csv'))
titles = total.columns

# ------------------- Descriptive table with all risks -----------------------
txt = '% --------------- written by risks.py from risks.csv ---------------\n'
txt += r'\begin{center}' + '\n' + r'\begin{longtable}{p{0.5cm}p{1.5cm}p{14cm}}' + '\n' + \
       r'    \caption{Technical risks.}\\' + '\n' + r'    \toprule' + '\n' + \
       r'    \multicolumn{1}{c}{\textbf{NB}}& \textbf{Risk ID} & \textbf{Description} \\' + '\n'

categories = total['category'].drop_duplicates()
for category in total['category'].drop_duplicates():
    txt += r'    \midrule' + '\n' + r'    \midrule' + '\n' + r'    \multicolumn{3}{c}{\textit{' + \
           str(category) + r'}} \\' + '\n' + r'    \midrule' + '\n'
    for idx, risk in total[total['category'] == str(category)].iterrows():
        txt += r'    \midrule' + '\n    ' + str(risk['number']) + ' & ' + str(risk['id']) + \
               r' & \textbf{' + str(risk['title']) + '} ' + str(risk['description']) + r' \newline \makecell{' + \
               '\n' + r'    \textit{Likelihood:} ' + str(risk['L1']) + ' |' + '\n' + \
               r'    \textit{Consequence:} ' + str(risk['C1']) + r'} \\' + '\n'

txt += r'    \bottomrule' + '\n' + r'    \label{tab:technical_risks}' + '\n' + r'\end{longtable}' + '\n' + \
       r'\end{center}' + '\n\n'

# ---------------------------- Construct pre-mitigation risk map ----------------------------
txt += r'\renewcommand{\arraystretch}{1.1} %stretches row height' + '\n' + r'\begin{table}[tb]' + '\n' + \
       r'\centering' + '\n' + r'\caption{\label{tab:risk_pre_mitigation} Technical risk matrix pre-mitigation.}' + \
       r'\resizebox{\textwidth}{!}' + '\n' + r'{\begin{tabular}{|p{2.2cm}|p{1.6cm}|p{1.9cm}|p{2.4cm}|p{2.7cm}|p{2.5cm}|}' + \
       r'\hline' + '\n' + r'& \textbf{Negligible} & \textbf{Minor} & \textbf{Moderate} & \textbf{Significant} & \textbf{Severe} \\ \hline' + '\n'

for likelihood in ['Very Likely', 'Likely', 'Possible', 'Unlikely', 'Very Unlikely']:
    txt += r'\textbf{' + likelihood + '}'
    for consequence in ['Negligible', 'Minor', 'Moderate', 'Significant', 'Severe']:
        txt += ' & ' + get_colour(likelihood, consequence)
        n = 0  # determine whether a comma is needed
        for idx, risk in total.iterrows():
            if risk['L1'] == likelihood and risk['C1'] == consequence:
                if n != 0:
                    txt += ', '
                txt += str(risk['number'])
                n += 1
    txt += r' \\ \hline' + '\n'

txt += r'\end{tabular}' + '\n' + r'}' + '\n' + r'\end{table}' + '\n' + \
       r'\renewcommand{\arraystretch}{1} %stretches row height' + '\n'

txt += r'\section{Mitigation Strategies}\label{sec:mitigation}'+ '\n' + \
       r'The risks in the top right corner of the risk map, with the colours orange or red, ' \
       r'are critical, and should be mitigated. These mitigation strategies, as well as their impact on the likelihood ' \
       r'and consequence, are described in \autoref{tab:risk_mitigation}. The risk map post-mitigation is shown in ' \
       r'\autoref{tab:risk_post_mitigation}.' + '\n'
# ----------------------------- Mitigation -----------------------------------
txt += r'\begin{center}' + '\n' + r'\begin{longtable}{p{0.5cm}p{15cm}}' + '\n' + \
       r'    \caption{Mitigation of technical risks.}\\' + '\n' + r'    \toprule' + '\n' + \
       r'    \textbf{Risk NB}& \textbf{Description \& Mitigation Plan} \\' + '\n'

mitigated_risks = total[total['mitigation'] == 'yes']
categories = mitigated_risks['category'].drop_duplicates()
for category in mitigated_risks['category'].drop_duplicates():
    for idx, risk in mitigated_risks[mitigated_risks['category'] == str(category)].iterrows():
        txt += r'    \midrule' + '\n    ' + str(risk['number']) + r' & \textbf{' + str(risk['title']) + r':} ' + str(risk['plan']) + r' \newline \makecell{' + \
               '\n' + r'    \textit{Likelihood:} ' + str(risk['L1']) + r' $\rightarrow$ ' + str(risk['L2']) + r' |' + '\n' + \
               r'    \textit{Consequence:} ' + str(risk['C1']) + r' $\rightarrow$ ' + str(risk['C2']) + r'} \\' + '\n'

txt += r'    \bottomrule' + '\n' + r'    \label{tab:risk_mitigation}' + '\n' + r'\end{longtable}' + '\n' + \
       r'\end{center}' + '\n\n'


# ---------------------------- Construct post-mitigation risk map ----------------------------
txt += r'\renewcommand{\arraystretch}{1.1} %stretches row height' + '\n' + r'\begin{table}[tb]' + '\n' + \
       r'\centering' + '\n' + r'\caption{\label{tab:risk_post_mitigation} Technical risk matrix post-mitigation.}' + \
       r'\resizebox{\textwidth}{!}' + '\n' + r'{\begin{tabular}{|p{2.2cm}|p{1.6cm}|p{1.9cm}|p{2.4cm}|p{2.7cm}|p{2.5cm}|}' + \
       r'\hline' + '\n' + r'& \textbf{Negligible} & \textbf{Minor} & \textbf{Moderate} & \textbf{Significant} & \textbf{Severe} \\ \hline' + '\n'

for likelihood in ['Very Likely', 'Likely', 'Possible', 'Unlikely', 'Very Unlikely']:
    txt += r'\textbf{' + likelihood + '}'
    for consequence in ['Negligible', 'Minor', 'Moderate', 'Significant', 'Severe']:
        txt += ' & ' + get_colour(likelihood, consequence)
        n = 0  # determine whether a comma is needed
        for idx, risk in total.iterrows():
            if risk['L2'] == likelihood and risk['C2'] == consequence:
                if n != 0:
                    txt += ', '
                txt += str(risk['number'])
                n += 1
    txt += r' \\ \hline' + '\n'

txt += r'\end{tabular}' + '\n' + r'}' + '\n' + r'\end{table}' + '\n' + \
       r'\renewcommand{\arraystretch}{1} %stretches row height' + '\n'

txt += '% ------------------------------------------------------------------\n'

# ------------------------------ Write to txt --------------------------------
file = open('risks.txt', 'w')
file.write(txt)
file.close()
print(txt)
