.. Dependency Plotting Additions for Tasks
   Mladen Ivkovic, Sep 2020


.. _task_adding_to_plotting_tool:
.. highlight:: python



Adding your Task to the Dependency Plotting Tool
================================================

How to create a task dependency graph is described in :ref:`Analysis_Tools`.
By default, it should pick up and plot the new task. However, you might want to
customize the plotting a bit, e.g. if you're introducing a new 'family' of tasks
or you'd like them plotted as a cluster instead of individually. In both cases, we
need to modify the ``tools/plot_task_dependencies.py`` script.



Colouring in the Task Nodes
---------------------------

First, the script needs to identify the task types with which it is working.
To do so, it will check the task names, which are generated following the scheme
``taskname_subtaskname``, where ``taskname`` is defined in ``taskID_names`` and
``subtasknbame`` is defined in ``subtaskID_names`` in ``task.c``. In 
``tools/plot_task_dependencies.py``, you'll have to write a function that recognizes your 
task by its name, like is done for example for gravity::

    def taskIsGravity(name):
        """
        Does the task concern the gravity?

        Parameters
        ----------

        name: str
            Task name
        """
        if "gpart" in name:
            return True
        if "grav" in name:
            return True
    return False

and add the check to the function ``writeTask()``::

    if taskIsGravity(name):
        txt += "color=red3,"

Feel free to pick out a `nice color <http://graphviz.org/doc/info/colors.html>`_ for it :)








Adding Clusters
---------------

In certain cases it makes sense to group some tasks together, for example the self 
and pair tasks when computing hydro densities, gradients, or forces. To do this, 
you'll need to modify the function ``task_get_group_name`` in ``src/task.c``. The group
is determined by the task subtype, e.g.

.. code-block:: c

    case task_subtype_grav:
      strcpy(cluster, "Gravity");
      break;

But since the task type itself is also passed to the function, you could use that
as well if you really really need to. And that's it!
