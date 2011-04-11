#ifndef _TRASHSTACK_H
#define _TRASHSTACK_H

struct trashstack {
    SLIST_HEAD(trashstack_head, trashstack_node) head;
};

struct trashstack_node {
    SLIST_ENTRY(trashstack_node) nodes;
};

static inline bool trashstack_init (struct trashstack *s)
{
    SLIST_INIT(&s->head);
    return s;
}

static inline void trashstack_deinit (struct trashstack *s)
{
    return;
}

static inline bool trashstack_empty (const struct trashstack *s)
{
    return SLIST_EMPTY(&s->head);
}

static inline ssize_t trashstack_size (const struct trashstack *s)
{
    const struct trashstack_node *node;
    ssize_t n = 0;
    
    SLIST_FOREACH(node, &s->head, nodes) {
        n++;
    }

    return n;
}

void * trashstack_pop (struct trashstack *s)
{
    void *node = NULL;
    
    if (!trashstack_empty(s)) {
        node = SLIST_FIRST(&s->head);
        SLIST_REMOVE_HEAD(&s->head, nodes);
    }

    return node;
}

void trashstack_push (struct trashstack *s, void *val)
{
    struct trashstack_node *node = val;
    SLIST_INSERT_HEAD(&s->head, node, nodes);
}


#endif /* _TRASHSTACK_H */
