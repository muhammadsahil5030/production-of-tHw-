import ROOT

#---------------- *Open the ROOT files and get the histogram from each file* ----------

fIn1 = ROOT.TFile.Open("/home/msahil/work/madgraph/MG5_aMC_v2_9_15/template_thwminus_signal1/Events/run_01_decayed_1/unweighted_events.root")

h1 = fIn1.Get("hM_wboson")

h1.SetLineWidth(2)

#--------------------------- *Scale the histograms by the cross sections* ------------
"""
xsec = [6.869e-00]
h1.Scale(xsec[0] / h1.Integral())
h2.Scale(xsec[1] / h2.Integral())
h3.Scale(xsec[2] / h3.Integral())
"""
#--------------------------- *Create a new canvas and draw the histograms* ------------

canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
h1.SetLineColor(ROOT.kRed)
h1.SetStats(0)

#h1.GetXaxis().SetRangeUser(0, 100)  # Set x-axis limits for hist1
#h1.GetYaxis().SetRangeUser(0,2600)  # Set y-axis limits for hist1

h1.Draw("E hist")

#-------------------------------- *Add the title and axis labels* ---------------------

h1.SetTitle("")
h1.GetXaxis().SetTitle("wboson mass [GeV]")
h1.GetYaxis().SetTitle("Events")

tex1 = ROOT.TLatex(0.153,0.92,"t#bar{t}t#bar{t}w- signal")
tex1.SetNDC()
tex1.SetTextAngle(0)
tex1.SetTextFont(42)
tex1.SetTextSize(0.04)
tex1.SetTextAlign(11)
tex1.Draw()

tex2 = ROOT.TLatex(0.730,0.92,"\sqrt{S} = {13TeV}")
tex2.SetNDC()
tex2.SetTextAngle(0)
tex2.SetTextFont(42)
tex2.SetTextAlign(11)
tex2.SetTextSize(0.04)
tex2.Draw()


#--------------------------------------- *Add a legend* -------------------------------

legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.85)
legend.AddEntry(h1, "w- mass", "l")
legend.SetBorderSize(0)
legend.Draw()

#-------------------------------------- *Show the canvas* -----------------------------
#canvas.SetLogy()
canvas.Draw()
canvas.SaveAs("wboson_mass.png")
input("press enter to exit....")

#-------------------------------------------- *End* -----------------------------------
