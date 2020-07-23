// Own header
#include "plot.h"

#include "TColor.h"

using namespace std;

//////////////////////////////////////////////////////////////////////
// globals

/// language in plots
bool gGerman = false;

/**********************************************************************
 * configuration functions
 */

char * strdup_new(const char * text)
{
  // duplicate string with help of new
  char * temp =  new char[strlen(text)+1];
  strcpy(temp, text);
  return temp;
}

/*********************************************************************
 * standard option settings
 */

void setopt(TStyle * style) 
{
  // style options
  style->SetStatColor(10);    // white background for statistics box
  style->SetPadColor(19);     // pads with grey background
  style->SetOptTitle(0);      // don't show histogram title
//   style->SetTitleColor(10);   // title background white
  style->SetFillColor(10);    // fill everything with a white background
  style->SetFillStyle(1001);  // solid fill style
  style->SetCanvasBorderMode(0); // no borders in canvases please!
  style->SetCanvasBorderSize(0); // no borders in canvases please!
  style->SetCanvasColor(10);  // white canvases please
  style->SetPadBorderMode(0); // no borders in Pads please!
  style->SetPadLeftMargin(0.15);
  style->SetPadBottomMargin(0.15);
  style->SetPadRightMargin(0.03);
  style->SetPadTopMargin(0.03);
  style->SetPadBorderSize(0); // no borders in Pads please!
  style->SetPadColor(10);     // white pads please
  style->SetOptStat(0);       // don't show statistics
  style->SetOptFit(1111);     // show fit parameters, probab., chi^2
  style->SetErrorX(0);        // no horizontal error
  style->SetTitleOffset(1.5, "Y"); // title offset
  style->SetTitleOffset(1.5, "X"); // title offset
  style->SetLabelSize(0.05, "Y"); // label size
  style->SetLabelSize(0.05, "X"); // label size
  style->SetTitleXSize(0.05);  // title size
  style->SetTitleYSize(0.05);  // title size
  style->SetFrameBorderMode(0);
  style->SetFrameBorderSize(0);
  style->SetMarkerStyle(kFullDotLarge);
  style->SetMarkerSize(1.);

  TColor mWhite(999, 1., 1., 1., "mWhite");

}

void setopt(TCanvas * canvas)
{
  // set default options for canvas
  canvas->SetTopMargin(0.03);
  canvas->SetRightMargin(0.05);
  canvas->SetBottomMargin(0.15);
  canvas->SetLeftMargin(0.15);
  canvas->SetFillColor(999);
  canvas->SetFillStyle(1001);
  canvas->SetBorderSize(0);
  canvas->SetBorderMode(0);
}

void setopt(TVirtualPad * pad)
{
  // set default options for canvas
  pad->SetTopMargin(0.05);
  pad->SetRightMargin(0.04);
  pad->SetBottomMargin(0.25);
  pad->SetLeftMargin(0.1);
  pad->SetFillColor(999);
  pad->SetFillStyle(1001);
  pad->SetBorderSize(0);
  pad->SetBorderMode(0);
}

void setopt(TH1 * histo)
{
  // set histo default options
  histo->SetLabelSize(0.06, "X"); // title offset
  histo->SetLabelSize(0.06, "Y"); // title offset
  histo->SetTitleSize(0.06, "X");  // title size
  histo->SetTitleSize(0.06, "Y");  // title size
  histo->SetTitleOffset(1.8, "Y"); // title offset
  histo->SetTitleOffset(1.1, "X"); // title offset
//    histo->GetXaxis()->SetTitleColor(kBlack);
//    histo->GetYaxis()->SetTitleColor(kBlack);
//    histo->GetXaxis()->CenterTitle();
}

void setoptlarge(TH1 * histo)
{
  // set histo default options
  histo->SetLabelSize(0.2, "X"); // title offset
  histo->SetLabelSize(0.2, "Y"); // title offset
  histo->SetTitleSize(0.2, "X");  // title size
  histo->SetTitleSize(0.2, "Y");  // title size
  histo->SetTitleOffset(0.2, "Y"); // title offset
  histo->SetTitleOffset(1., "X"); // title offset
//    histo->GetXaxis()->SetTitleColor(kBlack);
//    histo->GetYaxis()->SetTitleColor(kBlack);
//    histo->GetXaxis()->CenterTitle();
}

void setopt(TLegend * leg)
{
  // set legend default option
  leg->SetBorderSize(0);
  leg->SetFillColor(999);
}

void setopt(TGraph * gr)
{
  TH1F * histo = gr->GetHistogram();
  if (histo != 0)
    setopt(histo);
  gr->SetMarkerStyle(8);
  gr->SetMarkerSize(1.1);
}

const Text_t * GetUnit(const char * hname)
{
  // get unit of histogram name for labelling y-axis automatically
  static Text_t * unit = 0;

  if (unit != 0) {
    delete[] unit;
    unit = 0;
  }

  // find position where to separate the two strings
  Text_t pos2 = *strchr(hname, '@');
  Text_t * pos = &pos2;

  if (pos == 0) {
    // no unit given, return null pointer
    return 0;
  }

  // unit given, go over separator
  pos++;

  // and skip leading blanks
  while (isspace(*pos))
    pos++;

  unit = strdup_new(pos);

  return unit;
}

const char * GetXTitle(const TH1F * histo)
{
  // create new string containing only x title of histogram name
  static char * title = 0;

  if (title) {
    delete[] title;
    title = 0;
  }

  // find position where to separate the two strings
  const char * pos = strchr(histo->GetTitle(), ':');

  if (pos == 0) {
    pos = histo->GetTitle();
  }
  else {
    // y-title given, go over separator
    pos++;

    // and skip leading blanks
    while (isspace(*pos))
      pos++;
  }

  const char * upos = strchr(histo->GetTitle(), '@');
  if (upos != 0) {
    const char * unit = GetUnit(histo->GetTitle());
    // unit found
    Int_t len = upos-pos;
    title = new char[len+strlen(unit)+4];
    strncpy(title, pos, len);
    title[len] = '\0';
    strcat(title, " [");
    strcat(title, unit);
    strcat(title, "]");
  }
  else {
    // no unit
    title = strdup_new(pos);
  }
    
  return title;
}

const char * GetYTitle(const TH1F * histo)
{
  // create new string containing only y title of histogram name
  static char * title = 0;  // the y-title
  if (title) {
    delete[] title;
    title = 0;
  }
  // find position where to separate the two strings
  const char * pos = strchr(histo->GetTitle(), ':');

  // get unit
  const char * unit = GetUnit(histo->GetTitle());

  if (pos == 0) {
    // no y-title given
    Float_t binsize = histo->GetBinWidth(1);
    // unit given
    if (unit) {
      if (gGerman) {
	title = new char[26+strlen(unit)];
	sprintf(title, "Ereignisse / %4.2f %s", binsize, unit);
      }
      else {
	title = new char[25+strlen(unit)];
	sprintf(title, "Events / %4.2f %s", binsize, unit);
      }
    }
    else {
      if (gGerman) {
	title = new char[26];
	sprintf(title, "Ereignisse / %4.2f", binsize);
      }
      else {
	title = new char[25];
	sprintf(title, "Events / %4.2f", binsize);
      }
    }
  }
  else {
    const int len = strlen(histo->GetTitle())-strlen(pos);
    title = new char[len+1];
    strncpy(title, histo->GetTitle(), len);
    title[len] = '\0';
  }

  return title;
}

void settitle(TH1F * histo)
{
  histo->SetXTitle(GetXTitle(histo));
  histo->SetYTitle(GetYTitle(histo));
}
