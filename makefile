commit:
	git commit -am "commit"
	git push

commitall:
	git add -A
	$(MAKE) commit

pull:
	git pull
