if [ -d dist ] 
then
	rm -R dist
fi

mkdir dist
mkdir dist/src

cp -R -p /Library/Audio/Plug-Ins/Components/autotalent.component dist/autotalent.component
cp -R -p /Library/Audio/Plug-Ins/VST/autotalent.vst dist/autotalent.vst
cp -p ./* dist/src
cp -R -p ./autotalent.xcodeproj dist/src/autotalent.xcodeproj

ditto -c -k dist autotalent-mac-$(date +%d-%m-%Y).zip
rm -R dist
