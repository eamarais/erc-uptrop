Installation
==================

To download the latest stable version of the cloud-slicing uptrop-py code, type:

.. code-block:: console
  $ git clone https://github.com/eamarais/erc-uptrop.git

This command clones the source code to a folder called erc-uptrop

To download GEOS-Chem Classic source code into a folder named something other than :file:`erc-uptrop`, supply the name of the folder at the end of the :command:`git clone` command. For example:

.. code-block:: console
  $ git clone https://github.com/eamarais/erc-uptrop.git my-code-dir

This will download the source code into :file:`my-code-dir` instead of :file:`erc-uptrop`.

Once the :command:`git clone` process starts, you should see output similar to this:

.. code-block:: text

  Cloning into 'erc-uptrop'...
  remote: Enumerating objects: 1568, done.
  remote: Counting objects: 100% (485/485), done.
  remote: Compressing objects: 100% (277/277), done.
  remote: Total 1568 (delta 322), reused 297 (delta 195), pack-reused 1083
  Receiving objects: 100% (1568/1568), 1.25 MiB | 10.35 MiB/s, done.
  Resolving deltas: 100% (1086/1086), done.
  
When the :command:`git clone` process has finished, get a directory listing:

.. code-block:: console

   $ ls -CF erc-uptrop/*
   
and you will see output similar to this:

.. code-block:: text

   docs/  environment.yml  LICENSE  README.md  tests/  uptrop/

Navigate into the :file:`erc-uptrop` folder and confirm that the code is in the main branch:

.. code-block:: console

    $ cd erc-uptrop/
    $ git branch
    
You will see output similar to this:

.. code-block:: text

    * main

If you plan to modify the code, either to add new features or fix bugs, create a branch to store these changes. 

To do so, type:

.. code-block:: console

   $ git branch feature/my-git-updates
   $ git checkout feature/my-git-updates
   
Instead of :file:`feature/my-git-updates`, you may choose a name that reflects
the nature of your updates (e.g. :file:`feature/fix_bug`, :file:`feature/add_compound` etc.) 
If you now type:

.. code-block:: console

   $ git branch
   
You will see that we are checked out onto the branch that you just created.

.. code-block:: text

   * feature/my-git-updates
   main


