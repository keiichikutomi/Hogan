//----------------------------------------------------------------------------
#ifndef SDIMainH
#define SDIMainH
//----------------------------------------------------------------------------
#include <ComCtrls.hpp>
#include <ExtCtrls.hpp>
#include <Buttons.hpp>
#include <StdCtrls.hpp>
#include <Dialogs.hpp>
#include <Menus.hpp>
#include <Controls.hpp>
#include <Forms.hpp>
#include <Graphics.hpp>
#include <Classes.hpp>
#include <Windows.hpp>
#include <System.hpp>
#include <ActnList.hpp>
#include <ImgList.hpp>
#include <StdActns.hpp>
#include <ToolWin.hpp>
//----------------------------------------------------------------------------
class TSDIAppForm : public TForm
{
__published:
    TOpenDialog *OpenDialog;
    TSaveDialog *SaveDialog;
    TToolBar *ToolBar1;
    TToolButton *ToolButton9;
    TToolButton *ToolButton1;
    TToolButton *ToolButton2;
    TToolButton *ToolButton3;
    TToolButton *ToolButton4;
    TToolButton *ToolButton5;
    TToolButton *ToolButton6;
    TActionList *ActionList1;
    TAction *FileNew1;
    TAction *FileOpen1;
    TAction *FileSave1;
    TAction *FileSaveAs1;
    TAction *FileExit1;
    TEditCut *EditCut1;
    TEditCopy *EditCopy1;
    TEditPaste *EditPaste1;
    TAction *HelpAbout1;
    TStatusBar *StatusBar;
    TImageList *ImageList1;
    TMainMenu *MainMenu1;
    TMenuItem *File1;
    TMenuItem *FileNewItem;
    TMenuItem *FileOpenItem;
    TMenuItem *FileSaveItem;
    TMenuItem *FileSaveAsItem;
    TMenuItem *N1;
    TMenuItem *FileExitItem;
    TMenuItem *Edit1;
    TMenuItem *CutItem;
    TMenuItem *CopyItem;
    TMenuItem *PasteItem;
    TMenuItem *Help1;
    TMenuItem *HelpAboutItem;void __fastcall FileNew1Execute(TObject *Sender);
        void __fastcall FileOpen1Execute(TObject *Sender);
        void __fastcall FileSave1Execute(TObject *Sender);
        void __fastcall FileExit1Execute(TObject *Sender);
        void __fastcall HelpAbout1Execute(TObject *Sender);
private:
public:
    virtual __fastcall TSDIAppForm(TComponent *AOwner);
};
//----------------------------------------------------------------------------
extern TSDIAppForm *SDIAppForm;
//----------------------------------------------------------------------------
#endif
